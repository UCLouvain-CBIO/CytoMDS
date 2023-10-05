# CytoMDS - Copyright (C) <2023>
# <Université catholique de Louvain (UCLouvain), Belgique>
#
#   Description and complete License: see LICENSE file.
#
# This program (CytoMDS) is free software:
#   you can redistribute it and/or modify it under the terms of the GNU General
# Public License as published by the Free Software Foundation,
# either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details (<http://www.gnu.org/licenses/>).

#' @title Calculate Earth Mover's distance between two flowFrames
#'
#' @param ff1           a flowCore::flowFrame
#' @param ff2           a flowCore::flowFrame
#' @param channels      which channels (integer index(ices) or character(s)):
#' - if it is a character vector, 
#' it can refer to either the channel names, or the marker names
#' - if it is a numeric vector, 
#' it refers to the indexes of channels in `ff1`
#' - if NULL all scatter and fluorescent channels of `ff1`
#' will be selected
#' @param checkChannels  if `TRUE`, will explicitly check that
#' all provided channels are present in both flowFrames
#' @param binSize  size of equal bins to approximate 
#' the marginal distributions.
#' @param minRange minimum value taken 
#' when approximating the marginal distributions
#' @param maxRange maximum value taken 
#' when approximating the marginal distributions
#' @param returnAll If `TRUE`, distributions and marginal distribution
#' distances are returned as well. Default = `FALSE`.
#' @param nEventsLCM used to make the distance calculation a true metric. 
#' It is the lowest common multiplier of the number of events of each ff.
#' If it is provided as `NULL`, then it is calculated based on provided ffs.
#' If it is provided as an integer, then it is used as target (in the context
#' of calculating a matrix of pairwise distance between several ffs).
#'
#' @return the Earth Mover's distance between `ff1` and `ff2`,
#' which is calculated by summing up all EMD approximates for
#' the marginal distributions of each channel
#' @importFrom CytoPipeline areSignalCols
#' @export
#' 
#' @examples
#' 
#' library(CytoPipeline)
#' 
#' data(OMIP021Samples)
#' 
#' # estimate scale transformations 
#' # and transform the whole OMIP021Samples
#' 
#' transList <- estimateScaleTransforms(
#'     ff = OMIP021Samples[[1]],
#'     fluoMethod = "estimateLogicle",
#'     scatterMethod = "linearQuantile",
#'     scatterRefMarker = "BV785 - CD3")
#' 
#' OMIP021Trans <- CytoPipeline::applyScaleTransforms(
#'     OMIP021Samples, 
#'     transList)
#' 
#' # distance with itself (all channels at once)
#' # => should return 0
#' dist0 <- getEMDDist(
#'     ff1 = OMIP021Trans[[1]],
#'     ff2 = OMIP021Trans[[1]])
#' 
#' # returning only distance, 2 channels
#' dist1 <- getEMDDist(
#'     ff1 = OMIP021Trans[[1]], 
#'     ff2 = OMIP021Trans[[2]], 
#'     channels = c("FSC-A", "SSC-A"))
#' 
#' # using only one channel, passed by marker name
#' dist2 <- getEMDDist(ff1 = OMIP021Trans[[1]], 
#'                     ff2 = OMIP021Trans[[2]], 
#'                     channels = c("BV785 - CD3"))
#' 
#' # using only one channel, passed by index
#' dist3 <- getEMDDist(ff1 = OMIP021Trans[[1]], 
#'                     ff2 = OMIP021Trans[[2]], 
#'                     channels = 10)
#' 
#' dist2 == dist3
#' 
getEMDDist <- function(
        ff1, 
        ff2, 
        channels = NULL,
        checkChannels = TRUE,
        binSize = 0.05,
        minRange = -10,
        maxRange = 10,
        returnAll = FALSE,
        nEventsLCM = NULL) {
                       
    if (!inherits(ff1, "flowFrame") || !inherits(ff2, "flowFrame")) {
        stop("both flowFrame objects should inherit from flowCore::flowFrame")
    }
    
    ffList <- list(ff1, ff2)
    
    if (is.null(channels)) {
        channels <- flowCore::colnames(ff1)[areSignalCols(ff1)]
    } else if (is.numeric(channels)) {
        channels <- flowCore::colnames(ff1)[channels]
    } else {
        channels <- vapply(
            channels, 
            FUN = function(ch) {
                flowCore::getChannelMarker(ff1, ch)$name
            },
            FUN.VALUE = "c")
    }
    if (checkChannels) {
        # check that all channels are present in both flow frames
        vapply(X = ffList, FUN = function(ff) {
            wrongCh <- which(! channels %in% flowCore::colnames(ff))
            if (length(wrongCh) > 0) {
                stop(
                    "found some channels that are non existent in flowFrame ",
                    flowCore::identifier(ff), ":",
                    channels[wrongCh])
            }
            return(TRUE)
        }, FUN.VALUE = TRUE)  
    }
    
    # for performance
    ffList <- lapply(ffList, FUN = function(ff) ff[, channels])
    
    breaks <- seq(
        minRange, 
        maxRange, 
        by = binSize)
        
    
    breaks <- round(breaks,12)
    
    distr <- lapply(
        X = ffList,
        FUN = function(ff) {
            #browser()
            expr <- flowCore::exprs(ff)[, channels, drop = FALSE]
            
            # discretize all marginal distributions
            # check that the range correctly spans all events
            vapply(
                colnames(expr),
                exprMat = expr,
                FUN = function(colName, exprMat) {
                    #browser()
                    counts <- graphics::hist(
                        exprMat[, colName], 
                        breaks = c(-Inf, breaks, Inf),
                        plot = FALSE)$counts
                    counts <- counts[-c(1,length(counts))]
                    if (sum(counts) != nrow(exprMat)) {
                        warning(
                            "for flowFrame: [", 
                            flowCore::identifier(ff), "] : \n",
                            "provided [minRange, maxRange] does not ",
                            "span all events for channel ", colName,
                            ": count(events) = ", sum(counts), 
                            "; nEvents = ", nrow(exprMat))
                    }
                    return(counts)
                },
                FUN.VALUE = rep(0., length(breaks)-1)
                   
            )
        })
                    
    
    
    distances <- rep(0., length(channels))
    names(distances) <- channels
    
    nA <- flowCore::nrow(ff1)
    nB <- flowCore::nrow(ff2)
    
    ratioA <- 1
    ratioB <- 1
    if (is.null(nEventsLCM)) {
        nEventsLCM <-  pracma::Lcm(nA, nB)  
        ratioA <- nEventsLCM / nA
        ratioB <- nEventsLCM / nB
    } else {
        if (nEventsLCM - floor(nEventsLCM) > 1e-12)
            stop("provided [nEventsLCM] should be an integer")
        ratioA <- nEventsLCM / nA
        if (ratioA - floor(ratioA) > 1e-12)
            stop("provided [nEventsLCM] should be a multiple of nEvents(ff1)")
        
        ratioB <- nEventsLCM / nB
        if (ratioB - floor(ratioB) > 1e-12)
            stop("provided [nEventsLCM] should be a multiple of nEvents(ff2)")
    }
    
    for (ch in channels) {
        
        wA <- distr[[1]][, ch, drop=FALSE]
        wA <- wA * ratioA
        wB <- distr[[2]][, ch, drop=FALSE]
        wB <- wB * ratioB
        locations <- breaks[-1] - binSize/2
        # distances[ch] <- 
        #   emdist::emdw(A = locations,
        #                wA = wA,
        #                B = locations,
        #                wB = wB)
        distances[ch] <- 
            transport::wasserstein1d(
                a = locations,
                wa = wA,
                b = locations,
                wb = wB)
        
        # make sure distance is a PRECISE multiple 
        # of elementary mass transportation cost
        
        # if (distances[ch] < 1e-12) {
        #   distances[ch] <- 0
        # } else {
        #   elemCost <- binSize / nEventsLCM
        #   distances[ch] <- round(distances[ch]/elemCost) *  elemCost
        # }  
        
    }
    
    globalDist <- sum(distances)
    
    if (returnAll) {
        return(list(
            distr = distr, 
            distances = distances,
            globalDist = globalDist))
    }
    else {
        return(globalDist)
    }
}


#' @title Calculate all pairwise Earth Mover's distances 
#' between flowFrames of a flowSet
#' @param fs a flowCore::flowSet
#' @param channels which channels (integer index(ices) or character(s)):
#' - if it is a character vector, 
#' it can refer to either the channel names, or the marker names
#' - if it is a numeric vector, 
#' it refers to the indexes of channels in `fs`
#' - if NULL all scatter and fluorescent channels of `fs`
#' will be selected
#' @param verbose if `TRUE`, output a message 
#' after each single distance calculation
#' @param ... additional parameters passed to `getEMDDist()`
#' @return a distance matrix (full symmetric with 0. diagonal)
#' @importFrom CytoPipeline areSignalCols
#' @export
#' 
#' @examples
#' 
#' library(CytoPipeline)
#' 
#' data(OMIP021Samples)
#' 
#' # estimate scale transformations 
#' # and transform the whole OMIP021Samples
#' 
#' transList <- estimateScaleTransforms(
#'     ff = OMIP021Samples[[1]],
#'     fluoMethod = "estimateLogicle",
#'     scatterMethod = "linearQuantile",
#'     scatterRefMarker = "BV785 - CD3")
#' 
#' OMIP021Trans <- CytoPipeline::applyScaleTransforms(
#'     OMIP021Samples, 
#'     transList)
#'     
#' # calculate pairwise distances using only FSC-A & SSC-A channels
#' pwDist <- getPairwiseEMDDist(
#'     fs = OMIP021Trans,
#'     channels = c("FSC-A", "SSC-A"))
#' 
getPairwiseEMDDist <- function(
        fs,
        channels = NULL,
        verbose = FALSE,
        ...){
        
    if(!inherits(fs, "flowSet")) {
        stop("fs object should inherit from flowCore::flowSet")
    }
    
    # check channels
    if (is.null(channels)) {
        channels <- flowCore::colnames(fs)[areSignalCols(fs[[1]])]
    } else if (is.numeric(channels)) {
        channels <- flowCore::colnames(fs)[channels]
    } else {
        channels <- vapply(
            channels, 
            FUN = function(ch) {
                flowCore::getChannelMarker(fs[[1]], ch)$name
            },
            FUN.VALUE = "c")
    }
    
    # check that all channels are present in flowSet
    wrongCh <- which(! channels %in% flowCore::colnames(fs))
    if (length(wrongCh) > 0) {
        stop(
            "found some channels that are non existent in flowSet: ",
            channels[wrongCh])
    }
    
    # take only the channels of interest for the following
    #browser()
    # for performance
    fs <- fs[,channels]
    
    nFF <- length(fs)
    if (nFF < 1) stop("empty flowSet passed")
    
    
    # lowest common multiplier of the number of events of all flowframes
    #browser()
    nEventsLCM <- NULL
    # if (trueMetric) {
    #   
    #   fs <- flowCore::fsApply(fs,
    #                    FUN = function(ff){
    #                      nEvents <- flowCore::nrow(ff)
    #                      nEventsToAdd <- 2^ceiling(log2(nEvents)) - nEvents
    #                      if (nEventsToAdd > 0) {
    #                        withr::local_seed(0)
    #                        mySample <- sample(1:nEvents, size = nEventsToAdd)
    #                        flowCore::exprs(ff) <-
    #                          rbind(
    #                            flowCore::exprs(ff),
    #                            flowCore::exprs(ff)[mySample, , drop = FALSE]
    #                          )  
    #                      }
    #                      ff
    #                    })
    #   
    #   nCorrectedEvents <- flowCore::fsApply(fs, FUN = flowCore::nrow)
    #   
    #   nEventsLCM <- Reduce(f = pracma::Lcm, x = nCorrectedEvents, init = 1)
    #   #nEventsLCM <- max(nCorrectedEvents)
    # }
    
    #browser()
    pwDist <- diag(0., nrow = nFF)
    for (i in seq_len(nFF)) {
        j <- 1
        while (j < i) {
            # note channels are already checked :-)
            pwDist[i,j] <- getEMDDist(
                fs[[i]], 
                fs[[j]], 
                channels = channels,
                checkChannels = FALSE,
                nEventsLCM = nEventsLCM,
                ...)
                                      
            if (verbose) {
                message(
                    "i = ", i, 
                    "; j = ", j, 
                    "; dist = ", round(pwDist[i,j], 12))  
            }
            
            j <- j+1  
        }
    }
    
    # copy lower diagonal to upper diagonal
    pwDist <- pwDist + t(pwDist)
    
    return(pwDist)
}

#' @title Calculate a summary statistic of some channels of 
#' all flowFrames of a flowSet 
#' @param fs a flowCore::flowSet
#' @param channels which channels (integer index(ices) or character(s)):
#' - if it is a character vector, 
#' it can refer to either the channel names, or the marker names
#' - if it is a numeric vector, 
#' it refers to the indexes of channels in `fs`
#' - if NULL all scatter and fluorescent channels of `fs`
#' will be selected
#' @param verbose if `TRUE`, output a message 
#' after each single distance calculation
#' @param statFUNs a list (possibly of length one) of
#' functions to call to calculate the statistics, or a simple function
#' @param ... additional parameters passed to `getEMDDist()`
#' @return a matrix of which the columns are the channel statistics 
#' for all flowFrames of the flowSet. 
#' @importFrom CytoPipeline areSignalCols
#' @export
#' 
#' @examples
#' 
#' library(CytoPipeline)
#' 
#' data(OMIP021Samples)
#' 
#' # estimate scale transformations 
#' # and transform the whole OMIP021Samples
#' 
#' transList <- estimateScaleTransforms(
#'     ff = OMIP021Samples[[1]],
#'     fluoMethod = "estimateLogicle",
#'     scatterMethod = "linearQuantile",
#'     scatterRefMarker = "BV785 - CD3")
#' 
#' OMIP021Trans <- CytoPipeline::applyScaleTransforms(
#'     OMIP021Samples, 
#'     transList)
#'
#' channelsOrMarkers <- c("FSC-A", "SSC-A", "BV785 - CD3")
#' 
#' # calculate mean for each 4 selected channels, for each 2 samples
#' 
#' channelMeans <- getChannelsSummaryStat(
#'     OMIP021Trans,
#'     channels = channelsOrMarkers,
#'     statFUNs = mean)
#'     
#' # calculate median AND std deviation
#' # for each 4 selected channels, for each 2 samples
#' 
#' channelMedians <- getChannelsSummaryStat(
#'     OMIP021Trans,
#'     channels = channelsOrMarkers,
#'     statFUNs = list(stats::median, stats::sd))
#'    
getChannelsSummaryStat <- function(
        fs,
        channels = NULL,
        statFUNs = stats::median,
        verbose = FALSE,
        ...){
                                   
    if(!inherits(fs, "flowSet")) {
        stop("fs object should inherit from flowCore::flowSet")
    }
    
    nStats <- length(statFUNs)
    if (nStats < 1) {
        stop("At least one stat function should be provided for calculation")
    }
    
    if (nStats > 1) {
        statFUNs <- lapply(statFUNs, FUN = match.fun)
    } else {
        # also create a list to have it uniform single vs more functions cases
        statFUNs <- list(match.fun(statFUNs))
    }
    
    
    # check channels
    if (is.null(channels)) {
        channels <- flowCore::colnames(fs)[areSignalCols(fs[[1]])]
    } else if (is.numeric(channels)) {
        channels <- flowCore::colnames(fs)[channels]
    } else {
        channels <- vapply(
            channels, 
            FUN = function(ch) {
                flowCore::getChannelMarker(fs[[1]], ch)$name
            },
            FUN.VALUE = "c")
    }
    
    # chooses a display name for each channel:
    # markerName if not NA, otherwise channelName
    
    channelNames <- channels
    for (i in seq_along(channels)) {
        channelMarker <- flowCore::getChannelMarker(fs[[1]], channels[i])$desc
        if (!is.null(channelMarker) && !is.na(channelMarker)){
            channelNames[i] <- channelMarker
        }
    }
    
    # check that all channels are present in flowSet
    wrongCh <- which(! channels %in% flowCore::colnames(fs))
    if (length(wrongCh) > 0) {
        stop(
            "found some channels that are non existent in flowSet: ",
            channels[wrongCh])
    }
    
    # take only the channels of interest for the following
    #browser()
    # for performance
    fs <- fs[,channels]
    
    nFF <- length(fs)
    if (nFF < 1) stop("empty flowSet passed")
    
    nChannels <- length(channels)
    
    statList <- list()
    for (fu in seq_along(statFUNs)) {
        if (verbose) {
            message(
                "computing statistical function ", fu,
                " per channel...")
        }
        statList[[fu]] <- vapply(
            channels,
            FUN = function(ch){
                chRes <- flowCore::fsApply(
                    fs,
                    FUN = function(ff){
                        statFUNs[[fu]](flowCore::exprs(ff)[, ch], 
                                na.rm = TRUE)
                    })
            },
            FUN.VALUE = numeric(length = nFF))
        
        colnames(statList[[fu]]) <- channelNames
        rowNames <- flowCore::pData(fs)$name
        if (!is.null(rowNames)) {
            rownames(statList[[fu]]) <- rowNames
        }
    }
    
    # if only one stat function => unlist to return one single matrix
    if (nStats == 1){
        statList <- statList[[1]]
    }
    
    statList
}

#' @title metric MDS projection of sample
#' @description Multi-dimensional scaling projection of samples, 
#' using a distance matrix as an input. 
#' The MDS algorithm is not the classical MDS 
#' (`cmdscale` alike, aka Torgerson's algorithm), 
#' but is the SMACOF algorithm for metric distances that are not 
#' necessarily euclidean.  
#' After having obtained the projections on the `nDim` dimensions, 
#' we always apply svd decomposition to visualize as first axes the ones that 
#' contain the most variance of the projected dataset in `nDim` dimensions.   
#' Instead of being provided directly by the user, the `nDim` parameter can 
#' otherwise be found iteratively by finding the minimum `nDim` parameter that
#' allows the projection to reach a target pseudo RSquare.    
#' If this is the case, the `maxDim` parameter is used to avoid 
#' looking for too big projection spaces. 
#' @param pwDist (`nSamples` rows, `nSamples` columns), 
#' previously calculated pairwise distances between samples, 
#' must be provided as a full symmetric square matrix, with 0. diagonal
#' @param nDim number of dimensions of projection, as input to SMACOF algorithm
#' if not provided, will be found iteratively using `targetPseudoRSq`
#' @param seed seed to be set when launching SMACOF algorithm 
#' (e.g. when `init` is set to `"random"` but not only)
#' @param targetPseudoRSq target pseudo RSquare to be reached
#' (only used when `nDim` is set to NULL)
#' @param maxDim in case `nDim` is found iteratively, 
#' maximum number of dimensions the search procedure is allowed to explore
#' @param ... additional parameters passed to SMACOF algorithm
#'
#' @return a list with six elements:
#' - `$pwDist` the initial pair-wise distance (same as input)
#' - `$proj` the final configuration, i.e. the projected data matrix 
#' (`nSamples` rows, `nDim` columns) in `nDim` dimensions
#' - `$projDist` the distance matrix of projected data
#' - `stress` the global stress loss function final value 
#' obtained from the SMACOF algorithm
#' - `spp` the stress per point obtained from the SMACOF algorithm, i.e.
#' the contribution of each point to the stress loss function
#' - `$RSq` R squares, for each d, from 1 to `nDim`:
#' the (pseudo) R square when taking all dims from 1 to d.
#' - `$GoF` Goodness of fit, for each d, from 1 to `nDim`:
#' the goodness of fit indicator (b/w 0 and 1) when taking all dims from 1 to d.
#' Note pseudo R square and goodness of fit indicators are essentially the 
#' same indicator, only the definition of total sum of squares differ:
#' - for pseudo RSq: TSS is calculated using the mean pairwise distance 
#' as minimum
#' - for goodness of fit: TSS is calculated using 0 as minimum
#' 
#' @importFrom stats as.dist dist
#' @export
#' 
#' @examples
#' 
#' library(CytoPipeline)
#' 
#' data(OMIP021Samples)
#' 
#' # estimate scale transformations 
#' # and transform the whole OMIP021Samples
#' 
#' transList <- estimateScaleTransforms(
#'     ff = OMIP021Samples[[1]],
#'     fluoMethod = "estimateLogicle",
#'     scatterMethod = "linearQuantile",
#'     scatterRefMarker = "BV785 - CD3")
#' 
#' OMIP021Trans <- CytoPipeline::applyScaleTransforms(
#'     OMIP021Samples, 
#'     transList)
#'     
#' ffList <- flowCore::flowSet_to_list(OMIP021Trans)
#' 
#' # As there are only 2 samples in OMIP021Samples dataset,
#' # we create artificial samples that are random combinations of both samples
#' 
#' for(i in 3:5){
#'     ffList[[i]] <- 
#'         CytoPipeline::aggregateAndSample(
#'             OMIP021Trans,
#'             seed = 10*i,
#'             nTotalEvents = 5000)[,1:22]
#' }
#' 
#' fsNames <- c("Donor1", "Donor2", paste0("Agg",1:3))
#' names(ffList) <- fsNames
#' 
#' fsAll <- as(ffList,"flowSet")
#' 
#' # calculate all pairwise distances
#' 
#' pwDist <- getPairwiseEMDDist(fsAll, 
#'                              channels = c("FSC-A", "SSC-A"),
#'                              verbose = FALSE)
#' 
#' # compute Metric MDS object with explicit number of dimensions
#' mdsObj <- computeMetricMDS(pwDist, nDim = 2, seed = 0)
#' 
#' #' # compute Metric MDS object by reaching a target pseudo RSquare
#' mdsObj <- computeMetricMDS(pwDist, seed = 0, targetPseudoRSq = 0.999)
#' mdsObj$nDim # should be 3
#' 

computeMetricMDS <- function(
        pwDist,
        nDim = NULL,
        seed = NULL,
        targetPseudoRSq = 0.95,
        maxDim = 128,
        ...){
        
    #browser()
    
    dimensions <- dim(pwDist)
    if (length(dimensions) != 2) {
        stop("pwDist should be a square matrix")
    }
    if(dimensions[1] != dimensions[2]) {
        stop("pwDist should be a square matrix")
    }
    if (!is.numeric(pwDist)) {
        stop("pwDist should be numeric")
    }
    
    # handle case when nDim not provided
    # in that case iteratively find it by aiming at target pseudoRSquare
    
    if (is.null(nDim)) {
        nDimHigh <- 1
        currentRsq <- 0.
        while(currentRsq < targetPseudoRSq) {
            nDimHigh <- nDimHigh * 2
            if (nDimHigh > maxDim) {
                warning("maxDim (=", maxDim, ") reached without ",
                        "reaching target pseudo rsquare")
            }
            obj <- computeMetricMDS(pwDist,
                                    nDim = nDimHigh,
                                    seed = seed,
                                    ...)
            currentRsq <- obj$RSq[obj$nDim]
        }
        
        nDimLow <- 1
        
        # now [nDimLow, nDimUp] is such that RSq(nDimLow) < target 
        # and RSq(nDimHigh) is such that RSq(nDimHigh) >= target
        
        nDimMid <- nDimHigh # proper initialization ;-)
        while ((nDimHigh - nDimLow) > 1) {
            nDimMid <- floor((nDimLow + nDimHigh) / 2)
            obj <- computeMetricMDS(pwDist,
                                    nDim = nDimMid,
                                    seed = seed,
                                    ...)
            currentRsq <- obj$RSq[obj$nDim]
            if (currentRsq >= targetPseudoRSq) {
                nDimHigh <- nDimMid
            } else {
                nDimLow <- nDimMid
            }
        }
        # final recalculation of the last MDS was not done with the right nDim
        if (nDimMid != nDimHigh) {
            obj <- computeMetricMDS(pwDist,
                                    nDim = nDimHigh,
                                    seed = seed,
                                    ...)
        }
        
        return(obj)
    }
    
    # one-off case: nDim is provided by the user
    
    nSamples <- dimensions[1]
    if (nDim > nSamples-1) {
        stop("nDim should be at most (nSamples-1)")
    }
    res <- list()  
    
    if (!is.null(seed)) {
        withr::with_seed(
            seed,
            mdsRes <- smacof::smacofSym(
                delta = pwDist,
                ndim = nDim,
                principal = FALSE,
                ...)
        )
    } else {
        mdsRes <- smacof::smacofSym(
            delta = pwDist,
            ndim = nDim,
            principal = FALSE,
            ...)
    }
    proj <- mdsRes$conf
    
    # apply svd decomposition (principal component analysis) 
    # on the obtained projections
    
    proj_svd <- svd(proj)
    proj <- proj %*% proj_svd$v
    
    # store eigenvalues and pct of variance
    res$eigen <- proj_svd$d * proj_svd$d
    res$pctvar <- res$eigen/sum(res$eigen)
    
    delta <- as.dist(pwDist)
    N <- length(delta)
    scaleFactor <- sqrt(sum(delta^2, na.rm = TRUE)) / sqrt(N)
    proj <- proj * scaleFactor
    
    computeRSquares <- function(
        distances, 
        projections,
        asInLinearRegression = TRUE) {
        delta <- as.dist(distances)
        nDim <- ncol(projections)
        
        RSq <- rep(0., nDim)
        
        RSS <- vapply(
            seq_len(nDim), 
            FUN = function(d){
                projDist <- dist(projections[,seq_len(d)])
                sum((delta - projDist)^2)
            },
            FUN.VALUE = 0.)
        
        if (asInLinearRegression) {
            TSS <- sum((delta - mean(delta, na.rm = TRUE))^2)    
        } else {
            TSS <- sum(delta^2)    
        }
        
        
        
        RSq <- 1-RSS/TSS
        
        RSq
    }
    
    RSq <- computeRSquares(
        distances = pwDist,
        projections = proj, 
        asInLinearRegression = TRUE)
    GoF <- computeRSquares(
        distances = pwDist,
        projections = proj, 
        asInLinearRegression = FALSE)
    
    res$pwDist <- as.dist(pwDist)
    res$proj <- proj
    res$projDist <- dist(proj)
    res$stress <- mdsRes$stress
    res$spp <- mdsRes$spp
    res$RSq <- RSq
    res$GoF <- GoF
    res$nDim <- nDim
    res$mdsObj <- mdsRes
    
    class(res) <- "mdsRes"

    res
}

# Function to populate a mdsRes object with biplot information
# Note this information is specific to the couple of projectionAxes
# that will be used for the biplot
# returns a `mdsBiplot` list with 2 entries: 
# `$R2vec` and `$coefficients`
computeMetricMDSBiplot <- function(
        mdsObj,
        projectionAxes,
        extVariables) {
        
    X <- mdsObj$proj[,c(projectionAxes[1], projectionAxes[2])]
    p <- ncol(X)
    extVariables <- as.data.frame(extVariables)
    if (nrow(extVariables) != nrow(X)) {
        stop(
            "Number of rows in extVariables needs to match number of objects",
            "in configuration!")
    } 
    
    rownames(extVariables) <- rownames(mdsObj$proj)
    ext <- scale(extVariables, scale = TRUE)
    regfit <- lm(ext ~ -1 + X)
    mdsBiplot <- list()
    mdsBiplot$coefficients <- as.matrix(regfit$coefficients)
    regsum <- summary(regfit)
    if (ncol(ext) == 1) 
        
        R2vec <- regsum$r.squared
    else
        R2vec <- vapply(
            regsum,
            FUN = `[[`,
            FUN.VALUE = 0.,
            "r.squared")
    names(R2vec) <- colnames(ext)
    mdsBiplot$R2vec <- R2vec
    return(mdsBiplot)
}





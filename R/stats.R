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

# internal function calculating the unidimensional histograms of a flowFrame
.unidimHistograms <- function(
    ff, 
    breaks){
    
    expr <- flowCore::exprs(ff)#[, drop = FALSE]
    
    # discretize all marginal distributions
    # check that the range correctly spans all events
    distr <- vapply(
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
        FUN.VALUE = rep(0., length(breaks)-1))
    
    return(distr)
}

# internal function calculating the EMD distance from unidimensional histograms
.distFromUnidimHistograms <- function(breaks, distr1, distr2) {
    
    nBreaks <- length(breaks)
    if (nBreaks < 2) {
        stop("nBreaks should be at least 2")
    }
    channels <- colnames(distr1)
    nChannels <- length(channels)
    
    distances <- rep(0., length(channels))
    names(distances) <- channels
        
    nA <- sum(distr1[,1])
    nB <- sum(distr2[,1])
        
    ratioA <- 1
    ratioB <- 1
        
    nEventsLCM <-  pracma::Lcm(nA, nB)  
    ratioA <- nEventsLCM / nA
    ratioB <- nEventsLCM / nB
    
    #browser()
        
    for (ch in seq_along(channels)) {
        
        wA <- distr1[, ch, drop=FALSE]
        wA <- wA * ratioA
        wB <- distr2[, ch, drop=FALSE]
        wB <- wB * ratioB
        #locations <- breaks[-1] - binSize/2
        locations <- 0.5 * (breaks[-1] + breaks[-nBreaks])
        
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
    
    return(distances)
}

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
#' dist0 <- EMDDist(
#'     ff1 = OMIP021Trans[[1]],
#'     ff2 = OMIP021Trans[[1]])
#' 
#' # returning only distance, 2 channels
#' dist1 <- EMDDist(
#'     ff1 = OMIP021Trans[[1]], 
#'     ff2 = OMIP021Trans[[2]], 
#'     channels = c("FSC-A", "SSC-A"))
#' 
#' # using only one channel, passed by marker name
#' dist2 <- EMDDist(ff1 = OMIP021Trans[[1]], 
#'                     ff2 = OMIP021Trans[[2]], 
#'                     channels = c("BV785 - CD3"))
#' 
#' # using only one channel, passed by index
#' dist3 <- EMDDist(ff1 = OMIP021Trans[[1]], 
#'                     ff2 = OMIP021Trans[[2]], 
#'                     channels = 10)
#' 
#' dist2 == dist3
#' 
EMDDist <- function(
        ff1, 
        ff2, 
        channels = NULL,
        checkChannels = TRUE,
        binSize = 0.05,
        minRange = -10,
        maxRange = 10,
        returnAll = FALSE) {
                       
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
    ffList <- lapply(ffList, FUN = function(ff) ff[, channels, drop = FALSE])
    
    breaks <- seq(
        minRange,
        maxRange,
        by = binSize)
    
    breaks <- round(breaks,12)
    
    distrs <- lapply(
        X = ffList,
        # FUN = function(ff) {
        #     #browser()
        #     expr <- flowCore::exprs(ff)[, drop = FALSE]
        #     
        #     # discretize all marginal distributions
        #     # check that the range correctly spans all events
        #     vapply(
        #         colnames(expr),
        #         exprMat = expr,
        #         FUN = function(colName, exprMat) {
        #             #browser()
        #             counts <- graphics::hist(
        #                 exprMat[, colName], 
        #                 breaks = c(-Inf, breaks, Inf),
        #                 plot = FALSE)$counts
        #             counts <- counts[-c(1,length(counts))]
        #             if (sum(counts) != nrow(exprMat)) {
        #                 warning(
        #                     "for flowFrame: [", 
        #                     flowCore::identifier(ff), "] : \n",
        #                     "provided [minRange, maxRange] does not ",
        #                     "span all events for channel ", colName,
        #                     ": count(events) = ", sum(counts), 
        #                     "; nEvents = ", nrow(exprMat))
        #             }
        #             return(counts)
        #         },
        #         FUN.VALUE = rep(0., length(breaks)-1)
        #            
        #     )
        # }
        FUN = .unidimHistograms,
        breaks = breaks)
    
    #browser()
    
    
    # distances <- rep(0., length(channels))
    # names(distances) <- channels
    # 
    # nA <- flowCore::nrow(ff1)
    # nB <- flowCore::nrow(ff2)
    # 
    # ratioA <- 1
    # ratioB <- 1
    # 
    # nEventsLCM <-  pracma::Lcm(nA, nB)  
    # ratioA <- nEventsLCM / nA
    # ratioB <- nEventsLCM / nB
    # 
    # for (ch in channels) {
    #     
    #     wA <- distrs[[1]][, ch, drop=FALSE]
    #     wA <- wA * ratioA
    #     wB <- distrs[[2]][, ch, drop=FALSE]
    #     wB <- wB * ratioB
    #     locations <- breaks[-1] - binSize/2
    #     # distances[ch] <- 
    #     #   emdist::emdw(A = locations,
    #     #                wA = wA,
    #     #                B = locations,
    #     #                wB = wB)
    #     distances[ch] <- 
    #         transport::wasserstein1d(
    #             a = locations,
    #             wa = wA,
    #             b = locations,
    #             wb = wB)
    #     
    #     # make sure distance is a PRECISE multiple 
    #     # of elementary mass transportation cost
    #     
    #     # if (distances[ch] < 1e-12) {
    #     #   distances[ch] <- 0
    #     # } else {
    #     #   elemCost <- binSize / nEventsLCM
    #     #   distances[ch] <- round(distances[ch]/elemCost) *  elemCost
    #     # }  
    #     
    # }
    
    distances <- .distFromUnidimHistograms(
        breaks = breaks, 
        distr1 = distrs[[1]], 
        distr2 = distrs[[2]])
    
    globalDist <- sum(distances)
    
    if (returnAll) {
        return(list(
            distrs = distrs, 
            distances = distances,
            globalDist = globalDist))
    }
    else {
        return(globalDist)
    }
}

# internal function to generate blocks from a (rectangular piece of a) matrix
.generateBlocks2D <- function(
        rowRange, 
        colRange, 
        nRowBlock) {
    
    nRows <- rowRange[2] - rowRange[1] + 1
    nCols <- colRange[2] - colRange[1] + 1
    
    if (nRowBlock == 1) {
        blocks1DRows <- list()
        blocks1DRows[[1]] <- seq(rowRange[1], rowRange[2])
    } else {
        blocks1DRows <- split(seq(rowRange[1], rowRange[2]), 
                              cut(seq(rowRange[1], rowRange[2]), 
                                  nRowBlock, 
                                  labels = FALSE))
    }
    
    nColBlock <- ceiling(nRowBlock * nCols/nRows)
    
    if (nColBlock == 1) {
        blocks1DCols <- list()
        blocks1DCols[[1]] <- seq(colRange[1], colRange[2])
    } else {
        blocks1DCols <- split(seq(colRange[1], colRange[2]), 
                              cut(seq(colRange[1], colRange[2]), 
                                  nColBlock, 
                                  labels = FALSE))
    }
    
    blocks2D <- list()
    iBlocks2D <- 0
    for (j in seq_along(blocks1DRows)) {
        for (k in seq_along(blocks1DCols)) {
            if (blocks1DRows[[j]][1] < max(blocks1DCols[[k]])) {
                iBlocks2D <- iBlocks2D+1
                blocks2D[[iBlocks2D]] <- 
                    list(rowMin = blocks1DRows[[j]][1],
                         rowMax = max(blocks1DRows[[j]]),
                         colMin = blocks1DCols[[k]][1],
                         colMax = max(blocks1DCols[[k]]))
            }
        }
    }
    
    blocks2D
}   

# internal function to find a good 2D block size for a computation task,
# taking into account nb of available cores and memory size/core
.optimizeRowBlockNb <- function(
        rowRange, 
        colRange, 
        nCores = 1, 
        memSize = Inf) {
    
    nRows <- rowRange[2] - rowRange[1] + 1
    nCols <- colRange[2] - colRange[1] + 1
    
    # assumption: blocks are squares (apart from the borders) 
    # square side is unknown
    # => determine side of the square such that memSize is not exceeded
    # and nb of tasks is above nCores if possible.
    
    divisors <- seq(1, ceiling(sqrt(nRows)))
    candidateSides <- ceiling(nRows/divisors)
    candidateSides <- 
        union(candidateSides, seq(
            candidateSides[length(candidateSides)],
            1, by = -1))
    #candidateSides <- candidateSides[candidateSides<=sideMax]
    
    for (i in seq_along(candidateSides)) {
        nRowBlock <- ceiling(nRows/candidateSides[i])
        blocks2D <- .generateBlocks2D(
            rowRange, 
            colRange,
            nRowBlock = nRowBlock)
        if (length(blocks2D) == 0) {
            # special case obtained when there is only one element in
            # both rowRange and colRange
            nbInMem <- 0
        } else {
            nbInMem <- max(vapply(blocks2D,
                                  FUN = function(b){
                                      length(union(
                                          seq(b$rowMin, b$rowMax),
                                          seq(b$colMin, b$colMax)
                                      ))
                                  },
                                  FUN.VALUE = 1))    
        }
        
        if (nbInMem <= memSize && length(blocks2D) >= nCores) {
            #  memory threshold and suitable minimum nb of tasks reached
            break
        }
    }
    # note if the previous loop was not broken, this means
    # that finally each task will only compute one distance, and there are
    # not enough distances to compute to fill the nb of cores
    
    nRowBlock
}

# internal function for pairwise EMD distance calculation
# Note `memSize` has been hidden in the user interface function.
# It could be reactivated if we were to switch to distances not based on
# unidimensional histograms, in that case tight memory management will
# again be needed. 
#
# (3.)`memSize` should be an estimate of the number of flow frame (per core) 
# that can reside concurrently in memory. In that case the pairwise distance 
# calculation will be split by matrix blocs, and the block size 
# - hence the number of flow frames stored in memory concurrently - 
# will be adjusted according to the estimate.
# @param memSize specifies an estimate of the number of flowFrames that can
# live concurrently in the memory available to a single core 
# (in case BiocParallel is used). Note the provided value has to take into 
# account the type of BiocParallel infrastructure used (i.e. whether it uses 
# shared memory or not). 
.pairwiseEMDDist <- function(
        nSamples,
        rowRange = c(1, nSamples), 
        colRange = c(min(rowRange), nSamples),
        loadFlowFrameFUN,
        loadFlowFrameFUNArgs = NULL,
        channels = NULL,
        verbose = FALSE,
        useBiocParallel = FALSE,
        BPPARAM = BiocParallel::bpparam(),
        BPOPTIONS = BiocParallel::bpoptions(
            packages = c("flowCore")),
        memSize = Inf,
        binSize = 0.05,
        minRange = -10,
        maxRange = 10){
    
    # activate newer version of code
    # which works only when EMD distance is based on unidimensional
    # distribution distances
    unidimHistograms <- TRUE 
    
    # validate nSamples, rowSeq and colSeq
    if (!is.numeric(nSamples) || nSamples <= 0) {
        stop("nSamples should be an integer >= 0")
    }

    if (!is.numeric(rowRange) || length(rowRange) != 2) {
        stop("rowRange should be an integer vector of length 2")
    }
    if (rowRange[2] < rowRange[1]) {
        stop("rowRange should be ordered")
    }
    if (!all(seq(rowRange[1], rowRange[2]) %in% seq_len(nSamples))) {
        stop(
            "rowSeq elements should be included in the set of integers ",
            "<= nSamples")
    }

    if (!is.numeric(colRange) || length(colRange) != 2) {
        stop("colRange should be an integer vector of length 2")
    }
    if (colRange[2] < colRange[1]) {
        stop("rowRange should be ordered")
    }
    if (!all(seq(colRange[1], colRange[2]) %in% seq_len(nSamples))) {
        stop(
            "colSeq elements should be included in the set of integers ",
            "<= nSamples")
    }

    if (colRange[1] < rowRange[1]) {
        stop("rowRange and colRange should be such that the defined block is ",
             "from the upper triangular matrix (colRange[1] >= rowRange[1]")
    }
    
    if (!is.infinite(memSize)) {
        minMemSize <- ifelse(
            unidimHistograms,
            1,
            2)
        if (!is.numeric(memSize) || memSize < minMemSize) {
            stop("memSize should be an integer >= ", minMemSize)
        }
    }
    
    # from nb of available cores and possible memory restriction,
    # generate blocks to be run in one go
    if (useBiocParallel) {
        nAvailableCores <- BiocParallel::bpworkers(BPPARAM) 
    } else {
        nAvailableCores <- 1
    }
    
    #browser()
    
    # calculate the block2D task allocations
    nRowBlock <- .optimizeRowBlockNb(
        rowRange = rowRange,
        colRange = colRange,
        nCores = nAvailableCores,
        memSize = memSize)
    
    blocks2D <- .generateBlocks2D(
        rowRange = rowRange, 
        colRange = colRange, 
        nRowBlock = nRowBlock)
    
    if (unidimHistograms) {
        
        if (verbose){
            message("Pre-calculating all histograms...")
        }
        rowColSeqUnion <- union(
            seq(rowRange[1], rowRange[2]),
            seq(colRange[1], colRange[2]))
        
        # load all flowFrame in blocks and generate histograms
        if (nAvailableCores == 1){
            blocks1D <- list(rowColSeqUnion)
        } else {
            blocks1D <- split(rowColSeqUnion, 
                              cut(rowColSeqUnion, 
                              nAvailableCores, 
                              labels = FALSE)) 
        }
        
        breaks <- seq(
            minRange,
            maxRange,
            by = binSize)
        
        breaks <- round(breaks,12)
        
        loadFFAndCalcHistograms <- function(
            block, 
            loadFlowFrameFUN, 
            loadFlowFrameFUNArgs,
            channels,
            breaks,
            verbose) {
            
            #browser()
            
            distrs <- list()
            ind <- 0
            for (ffIndex in block){
                if (verbose) {
                    message("Loading file ", ffIndex, "...")
                }
                ind <- ind + 1
                ff <- do.call(
                    loadFlowFrameFUN,
                    args = c(list(ffIndex = ffIndex),
                             loadFlowFrameFUNArgs))
                invisible(gc()) # clean memory for ff stored in previous loop
                if(!inherits(ff, "flowFrame")) {
                    stop("object returned by loadFlowFrameFUN function ",
                         "should inherit from flowCore::flowFrame")
                }
                if (ind == 1) {
                    # check channels
                    if (is.null(channels)) {
                        channels <- 
                            flowCore::colnames(ff)[areSignalCols(ff)]
                    } else if (is.numeric(channels)) {
                        channels <- 
                            flowCore::colnames(ff)[channels]
                    } else {
                        channels <- vapply(
                            channels, 
                            FUN = function(ch) {
                                flowCore::getChannelMarker(ff, ch)$name
                            },
                            FUN.VALUE = "c")
                    }
                } else {
                    wrongCh <- which(! channels%in% flowCore::colnames(ff))
                    if (length(wrongCh) > 0) {
                        stop(
                            "found some channels that are non existent in flowframe #",
                            ffIndex, ": ",
                            channels[wrongCh])
                    } 
                }
                
                # take only the channels of interest for the following,
                # for performance
                ff <- ff[,channels]
                
                if (verbose) {
                    message("Calculating histogram for file ", ffIndex, "...")
                }
                distrs[[ind]] <- .unidimHistograms(
                    ff,
                    breaks = breaks
                )
            }
            distrs
        }
        
        if (useBiocParallel){
            distribBlockList <- BiocParallel::bplapply(
                blocks1D,
                FUN = loadFFAndCalcHistograms,
                BPPARAM = BPPARAM,
                BPOPTIONS = BPOPTIONS,
                loadFlowFrameFUN = loadFlowFrameFUN,
                loadFlowFrameFUNArgs = loadFlowFrameFUNArgs,
                channels = channels,
                breaks = breaks,
                verbose = verbose)
        } else {
            distribBlockList <- lapply(
                blocks1D,
                FUN = loadFFAndCalcHistograms,
                loadFlowFrameFUN = loadFlowFrameFUN,
                loadFlowFrameFUNArgs = loadFlowFrameFUNArgs,
                channels = channels,
                breaks = breaks, 
                verbose = verbose)
        }
        
        # reorganize multivariate distributions in a single list
        distrs <- list()
        ind <- 0
        for (b in seq_along(blocks1D)) {
            for(i in seq_along(blocks1D[[b]])) {
                ind <- ind+1    
                distrs[[ind]] <- distribBlockList[[b]][[i]]
            }
        }
        
        if (verbose){
            message("Calculating pairwise distances between histograms...")
        }
        
        handleOneBlockWithHistograms <- function(
            block,
            rowColSeqUnion, 
            breaks,
            distrs,
            verbose) {
            
            #browser()
            
            rowSeq <- seq(block$rowMin, block$rowMax)
            colSeq <- seq(block$colMin, block$colMax)
            nRows <- length(rowSeq)
            nCols <- length(colSeq)
            ffIndexes <- union(rowSeq, colSeq)
            
            pwDist <- matrix(rep(0., nRows * nCols),
                             nrow = nRows)
            
            for (i in seq_along(rowSeq)) {
                for (j in seq_along(colSeq)) {
                    if (colSeq[j] > rowSeq[i]) {
                        pwDist[i,j] <- sum(.distFromUnidimHistograms(
                            breaks = breaks,
                            distr1 = 
                                distrs[[which(rowColSeqUnion == rowSeq[i])]],
                            distr2 = 
                                distrs[[which(rowColSeqUnion == colSeq[j])]]))
                        if (verbose) {
                            message(
                                "i = ", rowSeq[i], 
                                "; j = ", colSeq[j], 
                                "; dist = ", round(pwDist[i,j], 12)) 
                        }
                    }
                }
            }
            
            # apply symmetry for block elements that are in the lower triangle
            for (i in seq_along(rowSeq)) {
                for (j in seq_along(colSeq)) {
                    if (colSeq[j] < rowSeq[i]) {
                        pwDist[i,j] <- pwDist[j,i]
                    }
                }
            }
            pwDist
        } 
        
        if (useBiocParallel) {
            pwDistByBlock <- BiocParallel::bplapply(
                blocks2D, 
                BPPARAM = BPPARAM,
                BPOPTIONS = BPOPTIONS,
                rowColSeqUnion = rowColSeqUnion,
                FUN = handleOneBlockWithHistograms,
                breaks = breaks,
                distrs = distrs,
                verbose = verbose)
        } else {
            pwDistByBlock <- lapply(
                blocks2D,
                rowColSeqUnion = rowColSeqUnion,
                FUN = handleOneBlockWithHistograms,
                breaks = breaks,
                distrs = distrs,
                verbose = verbose)
        }
    } else {
        handleOneBlock <- function(
        block,
        loadFlowFrameFUN,
        loadFlowFrameFUNArgs,
        channels,
        verbose) {
            
            rowSeq <- seq(block$rowMin, block$rowMax)
            colSeq <- seq(block$colMin, block$colMax)
            nRows <- length(rowSeq)
            nCols <- length(colSeq)
            ffIndexes <- union(rowSeq, colSeq)
            
            ffList <- list()
            
            ind <- 0
            for (ffIndex in ffIndexes){
                if (verbose) {
                    message("Loading file ", ffIndex, "...")
                }
                ind <- ind+1
                ff <- do.call(
                    loadFlowFrameFUN,
                    args = c(list(ffIndex = ffIndex),
                             loadFlowFrameFUNArgs))
                invisible(gc()) # clean memory for ff stored in previous loop
                if(!inherits(ff, "flowFrame")) {
                    stop("object returned by loadFlowFrameFUN function should inherit ",
                         "from flowCore::flowFrame")
                }
                if (ind == 1) {
                    # check channels
                    if (is.null(channels)) {
                        channels <- 
                            flowCore::colnames(ff)[areSignalCols(ff)]
                    } else if (is.numeric(channels)) {
                        channels <- 
                            flowCore::colnames(ff)[channels]
                    } else {
                        channels <- vapply(
                            channels, 
                            FUN = function(ch) {
                                flowCore::getChannelMarker(ff, ch)$name
                            },
                            FUN.VALUE = "c")
                    }
                } else {
                    wrongCh <- which(! channels%in% flowCore::colnames(ff))
                    if (length(wrongCh) > 0) {
                        stop(
                            "found some channels that are non existent in flowframe #",
                            ffIndex, ": ",
                            channels[wrongCh])
                    } 
                }
                
                # take only the channels of interest for the following,
                # for performance
                ffList[[ind]] <- ff[, channels]
            }
            
            pwDist <- matrix(rep(0., nRows * nCols),
                             nrow = nRows)
            
            for (i in seq_along(rowSeq)) {
                for (j in seq_along(colSeq)) {
                    if (colSeq[j] > rowSeq[i]) {
                        pwDist[i,j] <- 
                            EMDDist(
                                ff1 = ffList[[which(ffIndexes == rowSeq[i])]],
                                ff2 = ffList[[which(ffIndexes == colSeq[j])]],
                                channels = channels,
                                checkChannels = FALSE) 
                        if (verbose) {
                            message(
                                "i = ", rowSeq[i], 
                                "; j = ", colSeq[j], 
                                "; dist = ", round(pwDist[i,j], 12)) 
                        }
                    }
                }
            }
            
            # apply symmetry for block elements that are in the lower triangle
            for (i in seq_along(rowSeq)) {
                for (j in seq_along(colSeq)) {
                    if (colSeq[j] < rowSeq[i]) {
                        pwDist[i,j] <- pwDist[j,i]
                    }
                }
            }
            pwDist
        }
        
        if (useBiocParallel) {
            pwDistByBlock <- BiocParallel::bplapply(
                blocks2D, 
                BPPARAM = BPPARAM,
                BPOPTIONS = BPOPTIONS,
                FUN = handleOneBlock,
                loadFlowFrameFUN = loadFlowFrameFUN,
                loadFlowFrameFUNArgs = loadFlowFrameFUNArgs,
                channels = channels,
                verbose = verbose)
        } else {
            pwDistByBlock <- lapply(
                blocks2D,
                FUN = handleOneBlock,
                loadFlowFrameFUN = loadFlowFrameFUN,
                loadFlowFrameFUNArgs = loadFlowFrameFUNArgs,
                channels = channels,
                verbose = verbose)
        }
        
    }
    
    # sort out all block results to create one single matrix
    #browser()
    nRows <- rowRange[2] - rowRange[1] + 1
    nCols <- colRange[2] - colRange[1] + 1
    
    pwDist <- matrix(rep(0., nRows*nCols), nrow = nRows)
    for (b in seq_along(blocks2D)){
        block <- blocks2D[[b]]
        for (i in seq(block$rowMin, block$rowMax))
            for (j in seq(block$colMin, block$colMax))
                pwDist[i-rowRange[1]+1,j-colRange[1]+1] <- 
                    pwDistByBlock[[b]][i-block$rowMin+1,j-block$colMin+1]
    }
    # apply symmetry for lower triangular blocks
    rowSeq <- seq(rowRange[1], rowRange[2])
    colSeq <- seq(colRange[1], colRange[2])
    for (i in seq_along(rowSeq)) {
        for (j in seq_along(colSeq)){
            if (colSeq[j] < rowSeq[i]) {
                pwDist[i,j] <- pwDist[j,i]
            }
        }
    }
    
    rownames(pwDist) <- rowSeq
    colnames(pwDist) <- colSeq
    return(pwDist)
}

#' @title Pairwise Earth Mover's Distance calculation
#' @description Computation of all EMD between pairs of flowFrames belonging
#' to a flowSet.  
#' This method provides two different input modes:
#' - the user provides directly a flowSet loaded in memory (RAM).
#' - the user provides (1.) a number of samples `nSamples`; (2.) an ad-hoc 
#' function that takes as input an index between 1 and `nSamples`, and codes
#' the method to load the corresponding flowFrame in memory; 
#' Optional row and column ranges can be provided to limit the calculation
#' to a specific rectangle of the matrix. These i.e. can be specified as a way 
#' to split heavy calculations of large distance matrices 
#' on several computation nodes.  
#' 
#' @param x either a flowCore::flowSet, or the number of samples (integer >=1)
#' @param rowRange the range of rows of the distance matrix to be calculated
#' @param colRange the range of columns of the distance matrix to be calculated
#' @param loadFlowFrameFUN the function used to translate a flowFrame index
#' into a flowFrame. In other words, the function should code how to load a
#' specific flowFrame into memory. Important: the flow Frame index should be 
#' the first function argument and should be named `ffIndex`. 
#' @param loadFlowFrameFUNArgs (optional) a named list containing 
#' additional input parameters of `loadFlowFrameFUN()`
#' @param channels which channels (integer index(ices) or character(s)):
#' - if it is a character vector,
#' it can refer to either the channel names, or the marker names
#' - if it is a numeric vector,
#' it refers to the indexes of channels in `fs`
#' - if NULL all scatter and fluorescent channels of `fs` #' will be selected
#' @param verbose if `TRUE`, output a message
#' after each single distance calculation
#' @param useBiocParallel if `TRUE`, use `BiocParallel` for computation of the
#' pairwise distances in parallel - one (i,j) at a time.
#' Note the `BiocParallel` function used internally is `bplapply()`
#' @param BPPARAM if `useBiocParallel` is TRUE, sets the `BPPARAM` back-end to
#' be used for the computation. If not provided, will use the top back-end on
#' the `BiocParallel::registered()` stack.
#' @param BPOPTIONS if `useBiocParallel` is TRUE, sets the BPOPTIONS to be 
#' passed to `bplapply()` function.   
#' Note that if you use a `SnowParams` back-end, you need to specify all   
#' the packages that need to be loaded for the different CytoProcessingStep   
#' to work properly (visibility of functions). As a minimum,    
#' the `flowCore` package needs to be loaded.  
#' (hence the default `BPOPTIONS = bpoptions(packages = c("flowCore"))` )
#' @param binSize  size of equal bins to approximate 
#' the marginal distributions.
#' @param minRange minimum value taken 
#' when approximating the marginal distributions
#' @param maxRange maximum value taken 
#' when approximating the marginal distributions
#' @return a distance matrix of pairwise distances
#' (full symmetric with 0. diagonal)
#' @importFrom CytoPipeline areSignalCols
#' @importFrom methods as
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
#' pwDist <- pairwiseEMDDist(
#'     x = OMIP021Trans,
#'     channels = c("FSC-A", "SSC-A"))
#'
pairwiseEMDDist <- function(
        x,
        rowRange = c(1, nSamples), 
        colRange = c(min(rowRange), nSamples),
        loadFlowFrameFUN = NULL,
        loadFlowFrameFUNArgs = NULL,
        channels = NULL,
        verbose = FALSE,
        useBiocParallel = FALSE,
        BPPARAM = BiocParallel::bpparam(),
        BPOPTIONS = BiocParallel::bpoptions(
            packages = c("flowCore")),
        binSize = 0.05,
        minRange = -10,
        maxRange = 10
        ){
    
    nSamples <- NULL
    inMemory <- FALSE
    
    if(inherits(x, "flowSet")) {
        getFF <- function(ffIndex, fs) {
            return(fs[[ffIndex]])
        }
        nSamples <- length(x)
        pwDist <- .pairwiseEMDDist(
            nSamples = nSamples,
            rowRange = rowRange,
            colRange = colRange,
            loadFlowFrameFUN = getFF,
            loadFlowFrameFUNArgs = list(fs = x),
            channels = channels,
            verbose = verbose,
            useBiocParallel = useBiocParallel,
            BPPARAM = BPPARAM,
            BPOPTIONS = BPOPTIONS,
            binSize = binSize,
            minRange = minRange,
            maxRange = maxRange)
    } else {
        if (!is.numeric(x) || length(x) > 1) {
            stop("x should be either a flowCore::flowFrame ",
                 "or a numeric of length 1")
        }
        if (x < 1) {
            stop("x should be >= 1")
        }
        nSamples <- x
        
        if (is.null(loadFlowFrameFUN)) {
            stop("loadFlowFrameFUN should be provided ",
                 "when x is the number of samples")
        }
        
        if (useBiocParallel) {
            # blocks of maximum (approximately) 50*50 will be used
            memSize <- min(100, nSamples)
        } else {
            memSize <- Inf
        }
        
        pwDist <- .pairwiseEMDDist(
            nSamples = nSamples,
            rowRange = rowRange,
            colRange = colRange,
            loadFlowFrameFUN = loadFlowFrameFUN,
            loadFlowFrameFUNArgs = loadFlowFrameFUNArgs,
            channels = channels,
            verbose = verbose,
            useBiocParallel = useBiocParallel,
            BPPARAM = BPPARAM,
            BPOPTIONS = BPOPTIONS,
            memSize = memSize,
            binSize = binSize,
            minRange = minRange,
            maxRange = maxRange)
    }
    pwDist
}

# internal function calculating some by channel summary stats of a flowFrame
.calcFFStats <- function(
        ff, 
        channels,
        statFUNList,
        verbose){
    
    # chooses a display name for each channel:
    # markerName if not NA, otherwise channelName
    
    #browser()
    
    channelNames <- channels
    for (i in seq_along(channels)) {
        channelMarker <- flowCore::getChannelMarker(ff, channels[i])$desc
        if (!is.null(channelMarker) && !is.na(channelMarker)){
            channelNames[i] <- channelMarker
        }
    }

    nChannels <- length(channels)
    nStats <- length(statFUNList)
    
    statMatrix <- matrix(rep(0., nStats * nChannels), nrow = nStats)
    
    for (fu in seq_along(statFUNList)) {
        if (verbose) {
            message(
                "computing statistical function ", fu,
                " per channel...")
        }
        statMatrix[fu,] <- vapply(
            channels,
            FUN = function(ch){
                statFUNList[[fu]](
                flowCore::exprs(ff)[, ch, drop = FALSE],
                na.rm = TRUE)
            },
            FUN.VALUE = 0.)
    }
    
    colnames(statMatrix) <- channelNames
    rownames(statMatrix) <- names(statFUNList)
    
    
    statMatrix
}

# internal function for summary stats calculation
.channelSummaryStats <- function(
        nSamples,
        loadFlowFrameFUN,
        loadFlowFrameFUNArgs,
        channels = NULL,
        statFUNs = stats::median,
        verbose = FALSE,
        useBiocParallel = FALSE,
        BPPARAM = BiocParallel::bpparam(),
        BPOPTIONS = BiocParallel::bpoptions(
            packages = c("flowCore"))) {
    #browser()
    
    if (!is.numeric(nSamples) || nSamples < 1) {
        stop("nSamples should be a numeric >= 1")
    }
    
    nStats <- length(statFUNs)
    if (nStats < 1) {
        stop("At least one stat function should be provided for calculation")
    }
    
    if (is.list(statFUNs) || nStats>1) {
        statFUNList <- lapply(statFUNs, FUN = match.fun)
    } else {
        # also create a list to have it uniform single vs more functions cases
        statFUNList <- list(match.fun(statFUNs))
    }
    # set names to the list
    names(statFUNList) <- names(statFUNs)
    
    
    if (useBiocParallel) {
        nAvailableCores <- BiocParallel::bpworkers(BPPARAM) 
    } else {
        nAvailableCores <- 1
    }
    
    #browser()
    
    if (verbose){
        message("Loading flow frames and calculate stats...")
    }
    
    sequence <- seq_len(nSamples)
        
        
    if (nAvailableCores == 1){
        blocks1D <- list(sequence)
    } else {
        blocks1D <- split(sequence, 
                          cut(sequence, 
                              nAvailableCores, 
                              labels = FALSE)) 
    }
    
    loadFFAndCalcStats <- function(
        block, 
        loadFlowFrameFUN, 
        loadFlowFrameFUNArgs,
        channels,
        statFUNList,
        verbose) {
        
        #browser()
        
        statMatList <- list()
        ind <- 0
        for (ffIndex in block){
            if (verbose) {
                message("Loading file ", ffIndex, "...")
            }
            ind <- ind + 1
            ff <- do.call(
                loadFlowFrameFUN,
                args = c(list(ffIndex = ffIndex),
                         loadFlowFrameFUNArgs))
            invisible(gc()) # clean memory for ff stored in previous loop
            if(!inherits(ff, "flowFrame")) {
                stop("object returned by loadFlowFrameFUN function ",
                     "should inherit from flowCore::flowFrame")
            }
            if (ind == 1) {
                # check channels
                if (is.null(channels)) {
                    channels <- 
                        flowCore::colnames(ff)[areSignalCols(ff)]
                } else if (is.numeric(channels)) {
                    channels <- 
                        flowCore::colnames(ff)[channels]
                } else {
                    channels <- vapply(
                        channels, 
                        FUN = function(ch) {
                            flowCore::getChannelMarker(ff, ch)$name
                        },
                        FUN.VALUE = "c")
                }
            } else {
                wrongCh <- which(! channels%in% flowCore::colnames(ff))
                if (length(wrongCh) > 0) {
                    stop(
                        "found some channels that are non existent in flowframe #",
                        ffIndex, ": ",
                        channels[wrongCh])
                } 
            }
            
            # take only the channels of interest for the following,
            # for performance
            ff <- ff[,channels]
            
            if (verbose) {
                message("Calculating stats for file ", ffIndex, "...")
            }
            statMatList[[ind]] <- .calcFFStats(
                ff,
                channels = channels,
                statFUNList = statFUNList,
                verbose = verbose)
        } # end loop on ffIndex
        
        statMatList
        
    } # end function
    
    if (useBiocParallel){
        statMatBlockList <- BiocParallel::bplapply(
            blocks1D,
            FUN = loadFFAndCalcStats,
            BPPARAM = BPPARAM,
            BPOPTIONS = BPOPTIONS,
            loadFlowFrameFUN = loadFlowFrameFUN,
            loadFlowFrameFUNArgs = loadFlowFrameFUNArgs,
            channels = channels,
            statFUNList = statFUNList,
            verbose = verbose)
    } else {
        statMatBlockList <- lapply(
            blocks1D,
            FUN = loadFFAndCalcStats,
            loadFlowFrameFUN = loadFlowFrameFUN,
            loadFlowFrameFUNArgs = loadFlowFrameFUNArgs,
            channels = channels,
            statFUNList = statFUNList, 
            verbose = verbose)
    }
    
    #browser()
    
    # rearrange outputs
    
    nCh <- ncol(statMatBlockList[[1]][[1]])
    chStats <- list()
    for (s in seq_along(statFUNs)){
        chStats[[s]] <- matrix(
            data = rep(0., nSamples * nCh), nrow = nSamples)
        colnames(chStats[[s]]) <- colnames(statMatBlockList[[1]][[1]])
        ind <- 0
        for (b in seq_along(statMatBlockList)) {
            for (f in seq_along(statMatBlockList[[b]])) {
                ind <- ind + 1
                chStats[[s]][ind,] <- 
                    statMatBlockList[[b]][[f]][s,]
            }
        }
    }
    names(chStats) <- names(statFUNs)
    
    
    # if only one stat function => unlist to return one single matrix
    if (nStats == 1 && !is.list(statFUNs)){
        chStats <- chStats[[1]]
    }
    chStats
}

#' @title Summary statistics per channel computation
#' @description Computation of summary statistic for selected channels,
#' for all flowFrames of a flowSet.   
#' This method provides two different input modes:
#' - the user provides directly a flowSet loaded in memory (RAM).
#' - the user provides (1.) a number of samples `nSamples`; (2.) an ad-hoc 
#' function that takes as input an index between 1 and `nSamples`, and codes
#' the method to load the corresponding flowFrame in memory; 
#' Optional row and column ranges can be provided to limit the calculation
#' to a specific rectangle of the matrix. These i.e. can be specified as a way 
#' to split heavy calculations of large distance matrices 
#' on several computation nodes.  
#' @param x either a flowCore::flowSet, or the number of samples (integer >=1)
#' @param loadFlowFrameFUN the function used to translate a flowFrame index
#' into a flowFrame. In other words, the function should code how to load a
#' specific flowFrame into memory. Important: the flow Frame index should be 
#' the first function argument and should be named `ffIndex`. 
#' @param loadFlowFrameFUNArgs (optional) a named list containing 
#' additional input parameters of `loadFlowFrameFUN()`
#' @param channels which channels (integer index(ices) or character(s)):
#' - if it is a character vector, 
#' it can refer to either the channel names, or the marker names
#' - if it is a numeric vector, 
#' it refers to the indexes of channels in `fs`
#' - if NULL all scatter and fluorescent channels of `fs`
#' will be selected
#' @param statFUNs a list (possibly of length one) of
#' functions to call to calculate the statistics, or a simple function
#' This list can be named, in that case, these names will be transfered to the
#' returned value.
#' @param verbose if `TRUE`, output a message 
#' after each single distance calculation
#' @param useBiocParallel if `TRUE`, use `BiocParallel` for computation of the
#' pairwise distances in parallel - one (i,j) at a time.
#' Note the `BiocParallel` function used internally is `bplapply()`
#' @param BPPARAM if `useBiocParallel` is TRUE, sets the `BPPARAM` back-end to
#' be used for the computation. If not provided, will use the top back-end on
#' the `BiocParallel::registered()` stack.
#' @param BPOPTIONS if `useBiocParallel` is TRUE, sets the BPOPTIONS to be 
#' passed to `bplapply()` function.   
#' Note that if you use a `SnowParams` back-end, you need to specify all   
#' the packages that need to be loaded for the different CytoProcessingStep   
#' to work properly (visibility of functions). As a minimum,    
#' the `flowCore` package needs to be loaded.  
#' (hence the default `BPOPTIONS = bpoptions(packages = c("flowCore"))` )
#' @return a list of named statistic matrices. 
#' In each stat matrix, the columns are the channel statistics 
#' for all flowFrames of the flowSet.
#' Exception: if only one stat function (and not a list) is passed in
#' `statFUNs`, the return value is simplified to the stat matrix itself.
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
#' channelMeans <- channelSummaryStats(
#'     OMIP021Trans,
#'     channels = channelsOrMarkers,
#'     statFUNs = mean)
#'     
#' # calculate median AND std deviation
#' # for each 4 selected channels, for each 2 samples
#' 
#' channelMedians <- channelSummaryStats(
#'     OMIP021Trans,
#'     channels = channelsOrMarkers,
#'     statFUNs = list("median" = stats::median, 
#'                     "std.dev" = stats::sd))
#'  
channelSummaryStats <- function(
        x,
        loadFlowFrameFUN = NULL,
        loadFlowFrameFUNArgs = NULL,
        channels = NULL,
        statFUNs = stats::median,
        verbose = FALSE,
        useBiocParallel = FALSE,
        BPPARAM = BiocParallel::bpparam(),
        BPOPTIONS = BiocParallel::bpoptions(
            packages = c("flowCore"))){
    
    #browser()
    
    nSamples <- NULL
    inMemory <- FALSE
    
    if(inherits(x, "flowSet")) {
        getFF <- function(ffIndex, fs) {
            return(fs[[ffIndex]])
        }
        nSamples <- length(x)
        if (nSamples < 1) {
            stop("empty flowSet passed")
        }
        chStats <- .channelSummaryStats(
            nSamples = nSamples,
            loadFlowFrameFUN = getFF,
            loadFlowFrameFUNArgs = list(fs = x),
            channels = channels,
            statFUNs = statFUNs,
            verbose = verbose,
            useBiocParallel = useBiocParallel,
            BPPARAM = BPPARAM,
            BPOPTIONS = BPOPTIONS)
        
        # set flowframe names to returned matrix(ces)
        #browser()
        if (is.list(chStats)) {
            chStats <- lapply(chStats, FUN = function(e, fs){
                rownames(e) <- flowCore::sampleNames(fs)
                e
            }, fs = x)
        } else { # simple matrix
            rownames(chStats) <- flowCore::sampleNames(x)
        }
    } else {
        if (!is.numeric(x) || length(x) > 1) {
            stop("x should be either a flowCore::flowFrame ",
                 "or a numeric of length 1")
        }
        if (x < 1) {
            stop("x should be >= 1")
        }
        nSamples <- x
        
        if (is.null(loadFlowFrameFUN)) {
            stop("loadFlowFrameFUN should be provided ",
                 "when x is the number of samples")
        }
        
        chStats <- .channelSummaryStats(
            nSamples = nSamples,
            loadFlowFrameFUN = loadFlowFrameFUN,
            loadFlowFrameFUNArgs = loadFlowFrameFUNArgs,
            channels = channels,
            statFUNs = statFUNs,
            verbose = verbose,
            useBiocParallel = useBiocParallel,
            BPPARAM = BPPARAM,
            BPOPTIONS = BPOPTIONS)
    }
    chStats
    
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
#' flowCore::pData(fsAll)$type <- factor(c("real", "real", rep("synthetic", 3)))
#' flowCore::pData(fsAll)$grpId <- factor(c("D1", "D2", rep("Agg", 3)))
#' 
#' # calculate all pairwise distances
#' 
#' pwDist <- pairwiseEMDDist(fsAll, 
#'                              channels = c("FSC-A", "SSC-A"),
#'                              verbose = FALSE)
#' 
#' # compute Metric MDS object with explicit number of dimensions
#' mdsObj <- computeMetricMDS(pwDist, nDim = 4, seed = 0)
#' 
#' # mdsObj$nDim should be 4
#' 
#' #' # compute Metric MDS object by reaching a target pseudo RSquare
#' mdsObj2 <- computeMetricMDS(pwDist, seed = 0, targetPseudoRSq = 0.999)
#' 
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
# returns a `mdsBiplot` list with 3 entries: 
# `$R2vec`, `$coefficients` - the 2 latter obtained while regressing 
# the `extVariables` on the projectionAxes -, and `$correlations`, obtained
# by computing the Pearson correlation of the extVariables with the axes -
# each column being a vector of length nProjAxes corresponding to 
# the correlations of one ext variables with the projection axes.
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
    
    # calculate linear regressions
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
    
    # calculate Pearson correlations
    pearsonCorrelations <- t(stats::cor(ext, X))
    mdsBiplot$pearsonCorr <- pearsonCorrelations
    
    return(mdsBiplot)
}





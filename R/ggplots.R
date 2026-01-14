# CytoMDS - Copyright (C) <2023-2026>
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

#' @title Plot of channel intensity marginal densities
#' @description `ggplotMarginalDensities` uses ggplot2 
#' to draw plots of marginal densities of selected channels of a flowSet.
#' If the flowSet contains several flowFrames, events are concatenated 
#' together per group, or all together in the absence of groups.     
#' @param x a `flowCore::flowSet` (or a single `flowCore::flowFrame`)
#' @param sampleSubset (optional) a logical vector, of the same length as `x`, 
#' indicating which flow frames to keep in the plot. 
#' Typically it is obtained through the evaluation of a logical condition 
#' the rows of `phenoData(fs)`.   
#' @param channels (optional) - can be indices, or channel names, or markers.
#' @param pDataForColour (optional) which variable of `phenoData(fs)` 
#' will be used as colour aesthetic. Should be a character.
#' @param pDataForGroup (optional) which variable of `phenoData(fs)` 
#' will be used as group aesthetic. Should be a character.
#' A separate marginal density will be calculated for each group and 
#' overlaid on the same channel density plots.
#' @param nEventInSubsample how many event to take 
#' (per flowFrame of the flowSet) for marginal density approximation.
#' @param seed if not null, used in subsampling.
#' @param transList a `flowCore::transformList` that will be applied 
#' before plotting.
#' @export
#' @import ggplot2
#' @importFrom rlang .data
#' @return a ggplot object
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
#' 
#' # As there are only 2 samples in OMIP021Samples dataset,
#' # we create artificial samples that are random combinations of both samples
#' 
#' ffList <- c(
#'     flowCore::flowSet_to_list(OMIP021Trans),
#'     lapply(3:5,
#'            FUN = function(i) {
#'                aggregateAndSample(
#'                    OMIP021Trans,
#'                    seed = 10*i,
#'                    nTotalEvents = 5000)[,1:22]
#'            }))
#' 
#' fsNames <- c("Donor1", "Donor2", paste0("Agg",1:3))
#' names(ffList) <- fsNames
#' 
#' fsAll <- as(ffList,"flowSet")
#' 
#' flowCore::pData(fsAll)$grpId <- factor(c("D1", "D2", rep("Agg", 3)))
#' flowCore::pData(fsAll)$lbl <- paste0("S", 1:5)
#' 
#' # plot densities, all samples together
#' p <- ggplotMarginalDensities(fsAll)
#' 
#' # plot densities, per sample
#' p <- ggplotMarginalDensities(fsAll, pDataForGroup = "lbl")
#' 
#' # plot densities, per sample and coloured by group
#' p <- ggplotMarginalDensities(
#'     fsAll, 
#'     pDataForGroup = "lbl",
#'     pDataForColour = "grpId")
#' 
ggplotMarginalDensities <- function(
        x,
        sampleSubset,
        channels,
        pDataForColour, 
        pDataForGroup,
        nEventInSubsample = Inf,
        seed = NULL,
        transList) {
    
    if (inherits(x, "flowSet")){
        fs <- x
    } else if (inherits(x, "flowFrame")) {
        fs <- flowCore::flowSet(x)
    } else {
        stop("fs type not recognized, should be a flowSet or a flowFrame")
    }
    
    nSamples <- length(fs)
    if (nSamples == 0) {
        warning("empty flowSet passed")
    }
    
    if (missing(channels)) {
        ffAreSignalCols <- CytoPipeline::areSignalCols(fs[[1]])
        channels <- flowCore::colnames(fs[[1]])[ffAreSignalCols]
    }
    
    # sample subset
    if (missing(sampleSubset)) {
        sampleSubset <- rep_len(TRUE, nSamples)
    } else {
        if (!is.logical(sampleSubset) || length(sampleSubset) != nSamples) {
            stop("'sampleSubset' should be a logical vector of length = ",
                 "nb of samples")                            
        }
    }
    
    # transformation
    if (!missing(transList)) {
        if (!inherits(transList, "transformList")){
            stop("transList should be a tranformation list")
        }
        fs <- flowCore::transform(
            fs,
            transList)
    }
    
    
    if (nEventInSubsample < 1) 
        stop("nEventInSubsample should be strictly positive!")
    
    pData <- flowCore::pData(flowCore::phenoData(fs))
    pData <- pData[sampleSubset,, drop = FALSE]
    fs <- fs[sampleSubset]
    nSamples <- length(fs)
    
    # check pDataForColour and pDataForGroup
    if (!missing(pDataForColour)) {
        if (!is.character(pDataForColour)) {
            stop("'pDataForColour' should a character()")
        }
        if (!pDataForColour %in% flowCore::colnames(pData)) {
            stop("'pDataForColour' should be in phenoData columns")
        }
    } 
    if (!missing(pDataForGroup)) {
        if (!is.character(pDataForGroup)) {
            stop("'pDataForGroup' should a character()")
        }
        if (!pDataForGroup %in% flowCore::colnames(pData)) {
            stop("'pDataForGroup' should be in phenoData columns")
        }
    }
    
    nChannels <- length(channels)
    
    channelInfos <- lapply(channels,
                       FUN = function(ch, fr){
                           chmk <- flowCore::getChannelMarker(fr, ch)
                           if (is.null(chmk)) {
                               stop("channel ", ch, " not found in expr matrix")
                           }
                           name <- chmk$name
                           label <- name
                           if (!is.na(chmk$desc) && 
                               !toupper(chmk$desc) == "EMPTY"){
                               label <- chmk$desc
                           }
                           list(name = name, label = label)
                       },
                       fr = fs[[1]])
    
    channelNames <- vapply(channelInfos,
                           FUN = function(chi) chi$name,
                           FUN.VALUE = " ")
    channelLabels <- vapply(channelInfos,
                            FUN = function(chi) chi$label,
                            FUN.VALUE = " ")
    
    DFList <- mapply(
        seq_len(nSamples),
        FUN = function(i, fs, channels, channelLabels, nEventInSubsample, seed){
            nEvents <- flowCore::nrow(fs[[i]])
            chosenEvents <- 0
            if (nEventInSubsample < nEvents) {
                if (!is.null(seed)) {
                    withr::with_seed(
                        seed,
                        chosenEvents <- 
                            sample(nEvents, nEventInSubsample)
                    )
                } else {
                    chosenEvents <- sample(nEvents, nEventInSubsample)
                }
                nEvents <- nEventInSubsample
            } else {
                chosenEvents <- seq(nEvents)
            }
            DF <- data.frame(
                flowCore::exprs(fs[[i]])[chosenEvents, channels],
                pData[i,], 
                row.names = seq(chosenEvents))
            colnames(DF) <- c(channelLabels, colnames(pData))
            DF
        },
        MoreArgs = list(
            fs = fs,
            channels = channelNames,
            channelLabels = channelLabels,
            nEventInSubsample = nEventInSubsample,
            seed = seed),
        SIMPLIFY = FALSE)
    
    DF <- Reduce(DFList, f = rbind.data.frame)
    
    discardColourLegend <- FALSE
    if (missing(pDataForColour)) {
        DF$pDataForColour <- 1
        pDataForColour <- "pDataForColour"
        discardColourLegend <- TRUE
    }
    
    if (missing(pDataForGroup)) {
        DF$pDataForGroup <- 1
        pDataForGroup <- "pDataForGroup"
    }
    
    DFLong <- reshape2::melt(
        DF,
        measure.vars = channelLabels,
        value.name = "value",
        variable.name = "channel")
    
    p <- ggplot(
        DFLong, 
        aes(x = .data[["value"]],
            col = .data[[pDataForColour]],
            group = .data[[pDataForGroup]],
            y = after_stat(.data[["ndensity"]]))) + 
        facet_wrap(~channel, scales = "free_x") + 
        geom_density() + 
        ylab("normalized density")
    
    if (discardColourLegend){
        p <- p + guides(col = "none")
    }
    
    p 
} 

#' @title Plot of feature relative importance in distance
#' @description `ggplotDistFeatureImportance` uses ggplot2 
#' to provide a stacked bar plot of feature importance in a distance matrix.   
#' @param distObj a DistSum object.
#' @export
#' 
#' @seealso [pairwiseEMDDist]
#' 
#' @return a ggplot object
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
#' # As there are only 2 samples in OMIP021Samples dataset,
#' # we create artificial samples that are random combinations of both samples
#' 
#' ffList <- c(
#'     flowCore::flowSet_to_list(OMIP021Trans),
#'     lapply(3:5,
#'            FUN = function(i) {
#'                aggregateAndSample(
#'                    OMIP021Trans,
#'                    seed = 10*i,
#'                    nTotalEvents = 5000)[,1:22]
#'            }))
#' 
#' fsNames <- c("Donor1", "Donor2", paste0("Agg",1:3))
#' names(ffList) <- fsNames
#' 
#' fsAll <- as(ffList,"flowSet")
#' 
#' flowCore::pData(fsAll)$type <- factor(c("real", "real", rep("synthetic", 3)))
#' flowCore::pData(fsAll)$grpId <- factor(c("D1", "D2", rep("Agg", 3)))
#' 
#' # calculate all pairwise distances
#' 
#' pwDist <- pairwiseEMDDist(fsAll, 
#'                              channels = c("FSC-A", "SSC-A"),
#'                              verbose = FALSE)
#'                              
#' p <- ggplotDistFeatureImportance(pwDist)
#'                              
ggplotDistFeatureImportance <- function(distObj) {
    DF <- distByFeature(distObj)
    DF$featureName <- 
        factor(DF$featureName, 
               levels = DF$featureName[order(DF$percentage, decreasing = FALSE)])
    DF$dum <- factor(1)
    p <- ggplot(data = DF, aes(fill = .data[["featureName"]], 
                               x = .data[["dum"]], 
                               y = .data[["percentage"]])) + 
        geom_bar(stat = "identity", colour = "black") + 
        guides(fill = guide_legend(title = "feature")) + 
        theme(axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank()) + 
        labs(title = "Feature importance in distance matrix")
    
    p
}


#' @title Plot of Metric MDS object
#' @description `ggplotSampleMDS` uses ggplot2 
#' to provide plots of Metric MDS results.   
#' By default, a pseudo Rsquare projection quality indicator, 
#' and the number of dimensions of the MDS projection are provided in sub-title
#' @param mdsObj a MDS object, output of the `computeMetricMDS()` method.
#' @param pData (optional) a data.frame providing user input sample data. 
#' These can be design of experiment variables, phenotype data per sample,...
#' and will be used to highlight sample categories in the plot 
#' and/or for subsetting. 
#' @param sampleSubset (optional) a logical vector, of size `nrow(pData)`, 
#' which is by construction the nb of samples, indicating which samples to keep 
#' in the plot. Typically it is obtained through the evaluation of 
#' a logical condition on `pData` rows.   
#' @param projectionAxes which two axes should be plotted 
#' (should be a numeric vector of length 2)
#' @param biplot if TRUE, adds projection of external variables
#' @param extVariables are used to generate a biplot
#' these are the external variables that will be used in the biplot. 
#' They should be provided as a matrix with named columns corresponding to the 
#' variables. The number of rows should be the same as the number of samples.
#' The matrix might contain some NA's, in that case only complete rows will 
#' be used to calculate biplot arrows. 
#' @param biplotType type of biplot used:   
#' - if "correlation", projection of external variables will be according to 
#' Pearson correlations w.r.t. projection axes (arrow x & y coordinates)
#' - if "regression", a linear regression of external variables using the 2
#' projection axes as explanatory variables is performed, and the projection
#' of external variables will be according to regression coefficients
#' (arrow direction) and R square of regression (arrow size)
#' @param pDataForColour (optional) which `pData` variable
#' will be used as colour aesthetic. Should be a character.
#' @param pDataForShape (optional) which `pData` variable
#' will be used as shape aesthetic. Should be a character.
#' @param pDataForLabel (optional) which `pData` variable 
#' will be used as point labels in the plot. Should be a character.
#' If missing, point labels will be set equal to point names defined in 
#' MDS object (if not NULL, otherwise no labels will be set). 
#' @param pDataForAdditionalLabelling (optional) which `pData` variable(s)
#' will be added to the ggplot mapping, as to make them available for 
#' *plotly* tooltipping. Should be an array of character.
#' Note this works only if biplot=FALSE, as biplots contain circle and arrows 
#' that are currently not supported under `ggplotly`.
#' @param pointSizeReflectingStress if TRUE, size of points will appear 
#' proportional to stress by point, i.e. the bigger the sample point appears,
#' the less accurate its representation is 
#' (in terms of distances w.r.t. other points)
#' @param pointSize size of all points on the plots - only when   
#' `pointSizeReflectingStress` is FALSE. 
#' @param title title to give to the plot
#' @param displayPointLabels if TRUE, displays labels attached to points
#' (see `pDataForLabels` for the setting of the label values)
#' @param pointLabelSize size of point labels 
#' (default: 3.88 as in `geom_text()`)
#' @param repelPointLabels if TRUE, uses `ggrepel::geom_text_repel()` 
#' instead of `ggplot2::geom_text()`
#' (try to split the labels such that they do not overlap) for the points
#' @param displayArrowLabels if TRUE, displays arrows labels (only with biplot)
#' @param arrowLabelSize size of arrow labels 
#' (default: 3.88 as in `geom_text()`)  
#' @param repelArrowLabels if TRUE, uses `ggrepel::geom_text_repel()` 
#' instead of `ggplot2::geom_text()` for the arrows (only with biplot)
#' @param arrowThreshold (only with biplot), arrows will be made barely visible 
#' if their length is (in absolute value) less than this threshold.  
#' @param flipXAxis if TRUE, take the opposite of x values 
#' (provided as it might ease low dimensional projection comparisons)
#' @param flipYAxis if TRUE, take the opposite of y values 
#' (provided as it might ease low dimensional projection comparisons)
#' @param displayPseudoRSq if TRUE, display pseudo RSquare in subtitle, on top
#' of nb of dimensions
#' @param ... additional parameters passed to `ggrepel::geom_text_repel()` 
#' (if used)
#' @importFrom stats as.dist dist lm 
#' @importFrom rlang .data sym
#' 
#' @export
#' 
#' @seealso [ggplotSampleMDSWrapBiplots], [ggplotSampleMDSShepard], 
#' [computeMetricMDS]
#' 
#' @return a ggplot object
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
#' # As there are only 2 samples in OMIP021Samples dataset,
#' # we create artificial samples that are random combinations of both samples
#' 
#' ffList <- c(
#'     flowCore::flowSet_to_list(OMIP021Trans),
#'     lapply(3:5,
#'            FUN = function(i) {
#'                aggregateAndSample(
#'                    OMIP021Trans,
#'                    seed = 10*i,
#'                    nTotalEvents = 5000)[,1:22]
#'            }))
#' 
#' fsNames <- c("Donor1", "Donor2", paste0("Agg",1:3))
#' names(ffList) <- fsNames
#' 
#' fsAll <- as(ffList,"flowSet")
#' 
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
#' dim <- nDim(mdsObj) # should be 4
#' 
#' #' # compute Metric MDS object by reaching a target pseudo RSquare
#' mdsObj2 <- computeMetricMDS(pwDist, seed = 0, targetPseudoRSq = 0.999)
#' 
#' 
#' 
#' # plot mds projection on axes 1 and 2,
#' # use 'grpId' for colour, 'type' for shape, and no label 
#' 
#' p_12 <- ggplotSampleMDS(
#'     mdsObj = mdsObj,
#'     pData = flowCore::pData(fsAll),
#'     projectionAxes = c(1,2),
#'     pDataForColour = "grpId",
#'     pDataForShape = "type")
#' 
#' # plot mds projection on axes 3 and 4,
#' # use 'grpId' for colour, and 'name' as point label
#' 
#' p_34 <- ggplotSampleMDS(
#'     mdsObj = mdsObj,
#'     pData = flowCore::pData(fsAll),
#'     projectionAxes = c(3,4),
#'     pDataForColour = "grpId",
#'     pDataForLabel = "name")
#' 
#' # plot mds projection on axes 1 and 2,
#' # use 'group' for colour, 'type' for shape, and 'name' as point label
#' # have sample point size reflecting 'stress'
#' # i.e. quality of projection w.r.t. distances to other points
#' 
#' p12_Stress <- ggplotSampleMDS(
#'     mdsObj = mdsObj,
#'     pData = flowCore::pData(fsAll),
#'     projectionAxes = c(1,2),
#'     pDataForColour = "grpId",
#'     pDataForLabel = "name",
#'     pDataForShape = "type",
#'     pointSizeReflectingStress = TRUE)
#' 
#' # try to associate axes with median of each channel
#' # => use bi-plot
#' 
#' extVars <- channelSummaryStats(
#'     fsAll,
#'     channels = c("FSC-A", "SSC-A"),
#'     statFUNs = stats::median)
#' 
#' 
#' bp_12 <- ggplotSampleMDS(
#'     mdsObj = mdsObj,
#'     pData = flowCore::pData(fsAll),
#'     projectionAxes = c(1,2),
#'     biplot = TRUE,
#'     extVariables = extVars,
#'     pDataForColour = "grpId",
#'     pDataForShape = "type",
#'     seed = 0)
#' 
#' bp_34 <- ggplotSampleMDS(
#'     mdsObj = mdsObj,
#'     pData = flowCore::pData(fsAll),
#'     projectionAxes = c(3,4),
#'     biplot = TRUE,
#'     extVariables = extVars,
#'     pDataForColour = "grpId",
#'     pDataForLabel = "name",
#'     seed = 0)
#' 
ggplotSampleMDS <- function(
        mdsObj,
        pData,
        sampleSubset,
        projectionAxes = c(1,2),
        biplot = FALSE,
        biplotType = c("correlation", "regression"),
        extVariables,
        pDataForColour,
        pDataForShape,
        pDataForLabel,
        pDataForAdditionalLabelling,
        pointSize = 1.0,
        pointSizeReflectingStress = FALSE,
        title = "Multi Dimensional Scaling",
        displayPointLabels = TRUE,
        pointLabelSize = 3.88,
        repelPointLabels = TRUE,
        displayArrowLabels = TRUE,
        arrowLabelSize = 3.88,
        repelArrowLabels = FALSE,
        arrowThreshold = 0.8,
        flipXAxis = FALSE,
        flipYAxis = FALSE,
        displayPseudoRSq = TRUE,
        ...){
    
    if (!inherits(mdsObj, "MDS")) {
        stop("mdsObj should be a MDS object")
    }
    
    # biplot type (if any)
    biplotType <- match.arg(biplotType)
    
    nSamples <- nPoints(mdsObj)
    
    if (!missing(pData)) {
        if (nrow(pData) != nSamples) {
            stop("nb samples mismatch between mdsObj and pData")
        }
    }
    
    # sample subset
    if (missing(sampleSubset)) {
        sampleSubset <- rep_len(TRUE, nSamples)
    } else {
        if (!is.logical(sampleSubset) || length(sampleSubset) != nSamples) {
            stop("'sampleSubset' should be a logical vector of length = ",
                 "nb of samples")                            
        }
        # old code used for non standard evaluation (discarded)
        # e <- substitute(sampleSubset)
        # r <- eval(e, pData, parent.frame())
        # if (!is.logical(r)) 
        #     stop("'plotSubset' must be logical")
        # sampleSubset <- r & !is.na(r)    
    }
    
    if (!is.numeric(projectionAxes)) stop("projectionAxes should be numeric")
    if (!all(projectionAxes>=1)) stop("wrong values for projectionAxes")
    if (length(projectionAxes) != 2) stop("projectionAxes should have length 2")
    
    nDim <- nDim(mdsObj)
    
    if (nDim < max(projectionAxes)) {
        
        stop(
            "Nb of dimension of mds object (", 
            nDim, 
            ") is too low w.r.t. ",
            "requested projection axes")
    }
    
    if (biplot) {
        if (missing(extVariables)) {
            stop("biplot requires `extVariables`")
        }
        if (!is.numeric(extVariables)) {
            stop("`extVariables` should be a matrix")
        }
        dimensions <- dim(extVariables)
        if (length(dimensions) != 2) {
            stop("`extVariables` should be a matrix")
        }
        if (dimensions[1] != nSamples) {
            stop("`extVariables` should have nrow equal to nb of samples")
        }
    }
    
    # plot configuration
    
    # don't use notion of marginal RSquare anymore
    
    RSq <- RSq(mdsObj)
    GoF <- GoF(mdsObj)[nDim]
    
    explVar <- pctvar(mdsObj)
    
    proj <- projections(mdsObj)
    if (flipXAxis) {
        proj[, projectionAxes[1]] <- 
            - proj[, projectionAxes[1]]
    }
    if (flipYAxis) {
        proj[, projectionAxes[2]] <- 
            - proj[, projectionAxes[2]]
    }
    
    
    if (missing(pData)) {
        DF <- data.frame(
            x = proj[, projectionAxes[1]],
            y = proj[, projectionAxes[2]],
            stress = spp(mdsObj)
        )
    } else {
        DF <- pData
        DF$x <- proj[, projectionAxes[1]]
        DF$y <- proj[, projectionAxes[2]]
        DF$stress <- spp(mdsObj)
    }
    
    if (is.null(DF$sampleId)) {
        if (is.null(rownames(proj))) {
            DF$sampleId <- seq_len(nSamples)
        } else {
            DF$sampleId <- rownames(proj)
        }
    }
    
    if (missing(pDataForLabel)) {
        pDataForLabel <- "sampleId"
    }
    
    xlabel <- paste0("Coord. ", projectionAxes[1])
    
    xlabel <- paste0(
        xlabel, 
        " (% var. : ",
        round(100*explVar[projectionAxes[1]], 2),
        "%)")
    
    ylabel <- paste0("Coord. ", projectionAxes[2])
    
    ylabel <- paste0(
        ylabel, 
        " (% var. : ",
        round(100*explVar[projectionAxes[2]], 2),
        "%)")
    
    subtitle <- "("
    if (displayPseudoRSq) {
        subtitle <- paste0(
            subtitle,
            "Pseudo R2 = ", 
            round(RSq, 4),
            #"; Goodness of Fit = ", round(GoF, 4),
            "; ")
    } 
    subtitle <- paste0(
        subtitle, 
        "nDim = ",
        nDim,
        ")")
    
    # Character vectors for aesthetic keys and values
    mainAesKeys <- c("x", "y")
    mainAesValues <- c("x", "y")
    
    if (!missing(pDataForAdditionalLabelling) && !biplot) {
        if (!is.character(pDataForAdditionalLabelling)) {
            stop("pDataForAdditionalLabelling should be a character")
        }
        # remove from pDataForAdditionalLabelling :
        # 'x', 'y', pDataForColour, pDataForShape, pDataForLabel and 'stress'
        # to avoid possible duplicates
        toRemove <- c("x", "y")
        
        if (!missing(pDataForColour)) {
            toRemove <- union(toRemove, pDataForColour)
        }
        
        if (!missing(pDataForShape)) {
            toRemove <- union(toRemove, pDataForShape)
        }
        
        if (displayPointLabels && !missing(pDataForLabel)) {
            toRemove <- union(toRemove, pDataForLabel)
        }
        if (pointSizeReflectingStress) {
            toRemove <- union(toRemove, "stress")
        }
        
        pDataForAdditionalLabelling <- setdiff(pDataForAdditionalLabelling, 
                                               toRemove)
        
        if (length(pDataForAdditionalLabelling) > 0) {
            nAdditionalAes <- length(pDataForAdditionalLabelling) 
            mainAesKeys <- c(mainAesKeys, 
                             paste0("text", 1+seq_len(
                                 length(pDataForAdditionalLabelling))))
            mainAesValues <- c(mainAesValues, pDataForAdditionalLabelling)
        }
    }
    
    # dynamically create an aes mapping
    mainAesList <- lapply(mainAesValues, FUN = function(val) {
        rlang::sym(val)
    })
    names(mainAesList) <- mainAesKeys
    
    mainAesMapping <- do.call(aes, mainAesList)
    
    rangeMin <- min(proj)
    rangeMax <- max(proj)
    axesLimits <- c(rangeMin, rangeMax)
    #axesLimits <- c(min(-rangeMax, rangeMin), max(-rangeMin, rangeMax))
    
    p <- ggplot2::ggplot(
        data = DF[sampleSubset,],
        mapping = mainAesMapping) +
        ggplot2::labs(
            x = xlabel,
            y = ylabel,
            title = title,
            subtitle = subtitle) + 
        ggplot2::scale_x_continuous(limits = axesLimits) + 
        ggplot2::scale_y_continuous(limits = axesLimits)
    
    geomPointMapping <- NULL
    
    if (!missing(pDataForColour)) {
        if (!is.character(pDataForColour)) {
            stop("pDataForColour should be a character")
        }
        colourVar <- pDataForColour
        geomPointMapping <- c(
            geomPointMapping,
            ggplot2::aes(colour = .data[[colourVar]]))
    }
    if (!missing(pDataForShape)) {
        if (!is.character(pDataForShape)) {
            stop("pDataForShape should be a character")
        }
        shapeVar <- pDataForShape
        geomPointMapping <- c(
            geomPointMapping,
            ggplot2::aes(shape = .data[[shapeVar]]))
    }
    if (pointSizeReflectingStress) {
        geomPointMapping <- c(
            geomPointMapping,
            ggplot2::aes(size = .data[["stress"]]))
    }
    
    if (!is.null(geomPointMapping)) {
        # avoids error message: mapping should be created with `aes()`
        attr(geomPointMapping, "class") <- "uneval"
    }
    
    if (pointSizeReflectingStress) {
        p <- p + ggplot2::geom_point(mapping = geomPointMapping)    
    } else {
        p <- p + ggplot2::geom_point(
            mapping = geomPointMapping,
            size = pointSize)
    }
    
    if (displayPointLabels) {
        if (!is.character(pDataForLabel)) {
            stop("pDataForLabel should be a character")
        }
        labelVar <- pDataForLabel
        if (repelPointLabels) {
            p <- p + ggrepel::geom_text_repel(
                hjust=0.5, 
                vjust=1, 
                mapping = ggplot2::aes(label = .data[[labelVar]]),
                size = pointLabelSize,
                ...)
        } else {
            p <- p + ggplot2::geom_text(
                hjust=0.5, 
                vjust=1,
                mapping = ggplot2::aes(label = .data[[labelVar]]),
                size = pointLabelSize)
        }
    }
    
    # add biplot if specified
    if (biplot) {
        if (!is.numeric(arrowThreshold))
            stop("arrowThreshold should be numeric!")
        if (arrowThreshold < 0. || arrowThreshold > 1.) {
            stop("arrowThreshold should be between 0 and 1!")
        }
        
        mdsBiplot <- computeMetricMDSBiplot(
            mdsObj,
            projectionAxes = projectionAxes,
            extVariables = extVariables)
        
        radius <- 0.9*min(-axesLimits[1], axesLimits[2])
        lengthThreshold <- radius * arrowThreshold
        
        nExtVar <- ncol(mdsBiplot$coefficients)
        if (nExtVar > 0) {
            ret <- lapply(
                X = seq_len(nExtVar),
                FUN = function(j, mdsBiplot, radius, biplotType){
                    visible <- FALSE
                    valid <- FALSE
                    segmentName <- ""
                    segmentLength <- 0
                    segmentX <- 0
                    segmentY <- 0
                    if (biplotType == "regression"){
                        if (!is.na(mdsBiplot$R2vec[j])) {
                            valid <- TRUE
                            segmentLength <- radius * mdsBiplot$R2vec[j]
                            segmentAngle <- atan(
                                mdsBiplot$coefficients[2,j] / 
                                    mdsBiplot$coefficients[1,j])
                            if(mdsBiplot$coefficients[1,j] < 0){
                                segmentAngle <- segmentAngle + pi
                            }
                            segmentX <- segmentLength * cos(segmentAngle)
                            segmentY <- segmentLength * sin(segmentAngle)
                        }
                        
                    } else {
                        # biplotType is "correlation"
                        if (!is.na(mdsBiplot$pearsonCorr[1, j])) {
                            valid <- TRUE
                            segmentX <- radius * mdsBiplot$pearsonCorr[1, j]
                            segmentY <- radius * mdsBiplot$pearsonCorr[2, j]
                            segmentLength <- sqrt(segmentX^2 + segmentY^2)
                        }
                    }
                    
                    if (segmentLength >= lengthThreshold){
                        visible <- TRUE
                        segmentName <- colnames(mdsBiplot$coefficients)[j]
                    }
                    
                    
                    list(segmentX = segmentX,
                         segmentY = segmentY,
                         segmentName = segmentName,
                         visible = visible,
                         valid = valid)
                },
                mdsBiplot = mdsBiplot,
                radius = radius,
                biplotType = biplotType)
            
                
            # invert list hierarchy and simplify to array
            ret <- do.call(function(...) Map(list, ...), ret) 
            ret <- lapply(ret, FUN = simplify2array) 
            
            
            if (flipXAxis) {
                ret$segmentX <- - ret$segmentX  
            }
            if (flipYAxis) {
                ret$segmentY <- - ret$segmentY
            }
            
            segmentDF <- data.frame(
                segmentName = ret$segmentName,
                segmentXOrigin = rep(0., nExtVar),
                segmentYOrigin = rep(0., nExtVar),
                segmentX = ret$segmentX, 
                segmentY = ret$segmentY,
                visible = ret$visible)
            
            # discard non valid ext variables
            segmentDF <- segmentDF[ret$valid, ]
            
            p <- p + ggplot2::geom_segment(
                mapping = ggplot2::aes(
                    x = .data[["segmentXOrigin"]],
                    y = .data[["segmentYOrigin"]],
                    xend = .data[["segmentX"]],
                    yend = .data[["segmentY"]]),
                data = segmentDF[segmentDF$visible,],
                arrow = ggplot2::arrow(
                    length = ggplot2::unit(0.1, "inches"),
                    type = "closed"),
                linetype = 'dotted')
            
            newMapping <- ggplot2::aes(
                # x = (.data[["segmentXOrigin"]] + 
                #          .data[["segmentX"]]) / 2,
                # y = (.data[["segmentYOrigin"]] + 
                #          .data[["segmentY"]]) / 2
                x = .data[["segmentX"]],
                y = .data[["segmentY"]],
                hjust = "outward",
                vjust = "outward"
            )
            
            if (displayArrowLabels) {
                newMapping <- c(newMapping, 
                                ggplot2::aes(label = .data[["segmentName"]]))
            } else {
                # provide the segmentName for ggplotly
                newMapping <- c(newMapping,
                                ggplot2::aes(label = ""))
                #                      text1 = .data[["segmentName"]]))
            }
            
            # avoids error message: mapping should be created with `aes()`
            attr(newMapping, "class") <- "uneval"
            
            
            if (repelArrowLabels) {
                # discarding possible warning message: 
                # 'Ignoring unknown aesthetics: text1'
                p <- p + ggrepel::geom_text_repel(
                    mapping = newMapping,
                    data = segmentDF,
                    size = arrowLabelSize)
                
            } else {
                # discarding possible warning message: 
                # 'Ignoring unknown aesthetics: text1'
                p <- p + ggplot2::geom_text(
                    mapping = newMapping,
                    data = segmentDF,
                    size = arrowLabelSize)
            }
            
            # add a dashed centered circle to obtain the RSq=1 benchmark
            p <- p + ggforce::geom_circle(
                mapping = ggplot2::aes(
                    x0 = 0.,
                    y0 = 0.,
                    r = radius),
                linetype = 'dotted')
        } # if(nExtVar > 0)
    }
    
    p
}

#' @title Plot of Metric MDS object - Shepard diagram
#' @description `ggplotSampleMDSShepard` uses ggplot2 
#' to provide plot of Metric MDS results.  
#' Shepard diagram provides a scatter plot of :
#' - on the x axis, the high dimensional pairwise distances 
#' between each sample pairs  
#' - on the y axis, the corresponding pairwise distances in the obtained 
#' low dimensional projection
#' @param mdsObj a MDS object, output of the `computeMetricMDS()` method.
#' @param nDim (optional) number of dimensions to use when calculating   
#' Shepard's diagram and pseudoRSquare.  
#' If missing, it will be set equal to the number of projection dimensions  
#' as calculated in `mdsObj`
#' @param title title to give to the plot
#' @param pointSize point size in plot
#' @param lineWidth line width in plot 
#' @param displayPseudoRSq if TRUE, display pseudo RSquare in subtitle, on top
#' of nb of dimensions
#'
#' @importFrom rlang .data
#' @export
#' 
#' @seealso [ggplotSampleMDS], [computeMetricMDS]
#'
#' @return a ggplot object
#' 
#' @examples
#' 
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
#' ffList <- c(
#'     flowCore::flowSet_to_list(OMIP021Trans),
#'     lapply(3:5,
#'            FUN = function(i) {
#'                aggregateAndSample(
#'                    OMIP021Trans,
#'                    seed = 10*i,
#'                    nTotalEvents = 5000)[,1:22]
#'            }))
#' 
#' fsNames <- c("Donor1", "Donor2", paste0("Agg",1:3))
#' names(ffList) <- fsNames
#' 
#' fsAll <- as(ffList,"flowSet")
#' 
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
#' dim <- nDim(mdsObj) # should be 4
#' 
#' #' # compute Metric MDS object by reaching a target pseudo RSquare
#' mdsObj2 <- computeMetricMDS(pwDist, seed = 0, targetPseudoRSq = 0.999)
#' 
#' # Shepard diagrams 
#' 
#' p2D <- ggplotSampleMDSShepard(
#'     mdsObj,
#'     nDim = 2,
#'     pointSize = 1,
#'     title = "Shepard with 2 dimensions")
#' 
#' p3D <- ggplotSampleMDSShepard(
#'     mdsObj,
#'     nDim = 3,
#'     title = "Shepard with 3 dimensions") 
#'     #' 
#' pDefD <- ggplotSampleMDSShepard(
#'     mdsObj,
#'     title = "Shepard with default nb of dimensions") 
#'
ggplotSampleMDSShepard <- function(
        mdsObj,
        nDim,
        title = "Multi Dimensional Scaling - Shepard's diagram",
        pointSize = 0.5,
        lineWidth = 0.5,
        displayPseudoRSq = TRUE) {
    
    if (!inherits(mdsObj,"MDS")) {
        stop("mdsObj should be a MDS object")
    }
    
    if (missing(nDim)) {
        nDim <- CytoMDS::nDim(mdsObj)
    } else if (!is.numeric(nDim)) {
        stop("nDim should be numeric")
    } else if (nDim < 1) {
        stop("nDim should be >=1")
    } else if (nPoints(mdsObj) < nDim) {
        stop(
            "nDim too high compared to MDS object (nDim = ",
            nDim(mdsObj), ")")
    }
    
    RSq <- RSqVec(mdsObj)[nDim]
    GoF <- GoF(mdsObj)[nDim]
    
    subtitle <- "("
    if (displayPseudoRSq) {
        subtitle <- paste0(
            subtitle,
            "Pseudo R2 = ", 
            round(RSq, 4),
            #"; Goodness of Fit = ", round(GoF, 4),
            "; ")
    } 
    subtitle <- paste0(
        subtitle, 
        "nDim = ",
        nDim,
        ")")
    
    xlabel <- "HD distances"
    ylabel <- "Proj. distances"
    
    nSample <- nrow(pwDist(mdsObj))
    rowMatrix <- matrix(rep(seq_len(nSample), nSample), nrow = nSample)
    colMatrix <- matrix(
        rep(seq_len(nSample), nSample), 
        nrow = nSample, 
        byrow = TRUE)
    projDist <- as.vector(dist(projections(mdsObj)[,seq_len(nDim)]))
    HDDist <- as.vector(as.dist(pwDist(mdsObj)))
    rowNb <- as.vector(as.dist(rowMatrix))
    colNb <- as.vector(as.dist(colMatrix))
    
    DF <- data.frame(
        HDDist = HDDist,
        projDist = projDist,
        rowNb = rowNb,
        colNb = colNb)
    
    p <- ggplot2::ggplot(
        data = DF,
        mapping = ggplot2::aes(
            x = .data[["HDDist"]],
            y = .data[["projDist"]],
            text2 = .data[["rowNb"]],
            text3 = .data[["colNb"]])) +
        
        ggplot2::labs(
            x = xlabel,
            y = ylabel,
            title = title,
            subtitle = subtitle) + 
        
        ggplot2::geom_point(colour = "blue", size = pointSize) +
        ggplot2::geom_abline(
            intercept = 0.,
            slope = 1.,
            linetype = 'dotted',
            linewidth = lineWidth)
    p
}


#' @title SampleMDS biplot wrapping
#' @description `ggplotSampleMDSWrapBiplots` calls `ggplotSampleMDS` 
#' repeatly to generate biplots with different sets of external variables
#' and align them in a grid using the `patchwork` package, in a similar fashion 
#' as `ggplot2::facet_wrap()` does.
#' @param mdsObj a MDS object, output of the `computeMetricMDS()` method
#' @param extVariableList should be a named list of external variable matrices
#' Each element of the list should be a matrix with named columns 
#' corresponding to the variables. 
#' The number of rows should be the same as the number of samples. 
#' @param ncol passed to `patchwork::wrap_plots()`
#' @param nrow passed to `patchwork::wrap_plots()`
#' @param byrow passed to `patchwork::wrap_plots()`
#' @param displayLegend if FALSE, will de-active the legend display
#' @param ... additional parameters passed to `ggplotSampleMDS()` 
#' (if used)
#' 
#' @export
#' 
#' @seealso [ggplotSampleMDS], [ggplotSampleMDSShepard], [computeMetricMDS]
#' 
#' @return a ggplot object
#' 
#' @examples
#' 
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
#' # As there are only 2 samples in OMIP021Samples dataset,
#' # we create artificial samples that are random combinations of both samples
#' 
#' ffList <- c(
#'     flowCore::flowSet_to_list(OMIP021Trans),
#'     lapply(3:5,
#'            FUN = function(i) {
#'                aggregateAndSample(
#'                    OMIP021Trans,
#'                    seed = 10*i,
#'                    nTotalEvents = 5000)[,1:22]
#'            }))
#' 
#' fsNames <- c("Donor1", "Donor2", paste0("Agg",1:3))
#' names(ffList) <- fsNames
#' 
#' fsAll <- as(ffList,"flowSet")
#' 
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
#' dim <- nDim(mdsObj) # should be 4
#' 
#' #' # compute Metric MDS object by reaching a target pseudo RSquare
#' mdsObj2 <- computeMetricMDS(pwDist, seed = 0, targetPseudoRSq = 0.999)
#' 
#' # plot mds projection on axes 1 and 2,
#' # use 'group' for colour, 'type' for shape, and no label 
#' 
#' p_12 <- ggplotSampleMDS(
#'     mdsObj = mdsObj,
#'     pData = flowCore::pData(fsAll),
#'     projectionAxes = c(1,2),
#'     pDataForColour = "grpId",
#'     pDataForShape = "type")
#' 
#' # try to associate axes with median or std deviation of each channel
#' # => use bi-plots
#' 
#' extVarList <- channelSummaryStats(
#'     fsAll,
#'     channels = c("FSC-A", "SSC-A"),
#'     statFUNs = c("median" = stats::median, 
#'                  "std.dev" = stats::sd))
#' 
#' bpFull <- ggplotSampleMDSWrapBiplots(
#'     mdsObj = mdsObj,
#'     extVariableList = extVarList,
#'     pData = flowCore::pData(fsAll),
#'     projectionAxes = c(1,2),
#'     pDataForColour = "group",
#'     pDataForShape = "type",
#'     seed = 0)
#' 
ggplotSampleMDSWrapBiplots <- function(
        mdsObj,
        extVariableList,
        ncol = NULL,
        nrow = NULL,
        byrow = NULL,
        displayLegend = TRUE,
        ...){
    
    if (!is.list(extVariableList)) {
        stop("[extVariableList] should be a non zero length list!")
    }
    
    if (is.null(names(extVariableList))) {
        stop("[extVariableList] should be a named list!")
    }
    
    nPlots <- length(extVariableList)
    
    pList <- lapply(
        seq_len(nPlots),
        FUN = function(
            i, 
            extVariableList,
            mdsObj,
            displayLegend) {
            p <- ggplotSampleMDS(
                mdsObj = mdsObj,
                title = paste0("biplot with ",
                               names(extVariableList)[i]),
                biplot = TRUE,
                extVariables = extVariableList[[i]],
                ...) + ggplot2::labs(subtitle = NULL)
            if (!displayLegend) {
                p <- p + ggplot2:: theme(legend.position="none")
            }
            p
        },
        extVariableList = extVariableList,
        mdsObj = mdsObj,
        displayLegend = displayLegend
    )
        
    p <- patchwork::wrap_plots(
        pList, 
        ncol = ncol,
        nrow = nrow,
        byrow = byrow,
        guides = 'collect')
    p
}


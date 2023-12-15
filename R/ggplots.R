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


#' @title Plot of Metric MDS object
#' @description `ggplotSampleMDS` uses ggplot2 
#' to provide plots of Metric MDS results.   
#' By default, a pseudo Rsquare projection quality indicator, 
#' and the number of dimensions of the MDS projection are provided in sub-title
#' @param mdsObj a MDS object calculated by the SMACOF algorithm using
#' the computeMetricMDS() function
#' @param pData a data.frame providing user input sample data. 
#' These can be design of experiment variables, phenotype data per sample,...
#' and will be used to highlight sample categories in the plot. 
#' @param projectionAxes which two axes should be plotted 
#' (should be a numeric vector of length 2)
#' @param biplot if TRUE, adds projection of external variables
#' @param extVariables are used to generate a biplot
#' these are the external variables that will be used in the biplot. 
#' They should be provided as a matrix with named columns corresponding to the 
#' variables. The number of rows should be the same as the number of samples. 
#' @param biplotType type of biplot used:   
#' - if "correlation", projection of external variables will be according to 
#' Pearson correlations w.r.t. projection axes (arrow x & y coordinates)
#' - if "regression", a linear regression of external variables using the 2
#' projection axes as explanatory variables is performed, and the projection
#' of external variables will be according to regression coefficients
#' (arrow direction) and R square of regression (arrow size)
#' @param pDataForColour if not NULL, which `pData` variable
#' will be used as colour aesthetic. Should be a character.
#' @param pDataForShape if not NULL, which `pData` variable
#' will be used as shape aesthetic. Should be a character.
#' @param pDataForLabel if not NULL, which `pData` variable 
#' will be used as point labels in the plot. Should be a character.
#' @param pDataForAdditionalLabelling if not NULL, which `pData` variable(s)
#' will be add to the ggplot mapping, as to make them available for 
#' *plotly* tooltipping. Should be an array of character of maximum length 3.
#' Note this works only if biplot=FALSE, as biplots contain circle and arrows 
#' that are currently not supported under `ggplotly`.
#' @param sizeReflectingStress if TRUE, size of points will appear 
#' proportional to stress by point, i.e. the bigger the sample point appears,
#' the less accurate its representation is 
#' (in terms of distances w.r.t. other points)
#' @param title title to give to the plot
#' @param repelPointsLabels if TRUE, uses `ggrepel::geom_text_repel()` 
#' instead of `ggplot2::geom_text()`
#' (try to split the labels such that they do not overlap) for the points
#' @param repelArrowLabels if TRUE, uses `ggrepel::geom_text_repel()` 
#' instead of `ggplot2::geom_text()` for the arrows (only with biplot)
#' @param displayArrowLabels if TRUE, displays arrows labels (only with biplot)
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
#' @importFrom rlang .data
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
#' # prepare data, build MDS object
#' example("computeMetricMDS")
#' 
#' # plot mds projection on axes 1 and 2,
#' # use 'grpId' for colour, 'type' for shape, and no label 
#' 
#' p_12 <- ggplotSampleMDS(
#'     mdsObj = mdsObj,
#'     pData = flowCore::pData(fsAll),
#'     projectionAxes = c(1,2),
#'     pDataForColour = "grpId",
#'     pDataForLabel = NULL,
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
#'     sizeReflectingStress = TRUE)
#' 
#' # try to associate axes with median of each channel
#' # => use bi-plot
#' 
#' extVars <- getChannelsSummaryStat(
#'     fsAll,
#'     channels = c("FSC-A", "SSC-A"),
#'     statsFUN = stats::median)
#' 
#' 
#' bp_12 <- ggplotSampleMDS(
#'     mdsObj = mdsObj,
#'     pData = flowCore::pData(fsAll),
#'     projectionAxes = c(1,2),
#'     biplot = TRUE,
#'     extVariables = extVars,
#'     pDataForColour = "grpId",
#'     pDataForLabel = NULL,
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
        pData = data.frame(
            sampleId = seq_len(nrow(mdsObj$proj))),
        projectionAxes = c(1,2),
        biplot = FALSE,
        biplotType = c("correlation", "regression"),
        extVariables = NULL,
        pDataForColour = NULL,
        pDataForShape = NULL,
        pDataForLabel = "name",
        pDataForAdditionalLabelling = NULL,
        sizeReflectingStress = FALSE,
        title = "Multi Dimensional Scaling",
        repelPointsLabels = TRUE,
        displayArrowLabels = TRUE,
        arrowThreshold = 0.8,
        repelArrowLabels = FALSE,
        flipXAxis = FALSE,
        flipYAxis = FALSE,
        displayPseudoRSq = TRUE,
        ...){
    
    #browser()
    
    if (!inherits(mdsObj, "mdsRes")) {
        stop("mdsObj should be a 'mdsRes' object")
    }
    
    biplotType <- match.arg(biplotType)
    
    nSamples <- nrow(pData)
    
    if (nrow(mdsObj$proj) != nSamples) {
        stop("nb samples mismatch between mdsObj and pData")
    }
    
    if (!is.numeric(projectionAxes)) stop("projectionAxes should be numeric")
    if (!all(projectionAxes>=1)) stop("wrong values for projectionAxes")
    if (length(projectionAxes) != 2) stop("projectionAxes should have length 2")
    
    nDim <- ncol(mdsObj$proj)
    
    if (nDim < max(projectionAxes)) {
        
        stop(
            "Nb of dimension of mds object (", 
            nDim, 
            ") is too low w.r.t. ",
            "requested projection axes")
    }
    
    nSamples <- nrow(mdsObj$proj)
    
    if (biplot) {
        if (is.null(extVariables)) {
            stop("biplot requires non null `extVariables`")
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
    
    #browser()
    
    RSq <- mdsObj$RSq[nDim]
    GoF <- mdsObj$GoF[nDim]
    
    # margRSq <- rep(0., nDim)
    # 
    # for(j in 2:nDim) {
    #     margRSq[j] <- RSq[j] - RSq[j-1]
    # }
    # margRSq[1] <- RSq[1]
    
    explVar <- mdsObj$pctvar
    
    proj <- mdsObj$proj
    if (flipXAxis) {
        proj[, projectionAxes[1]] <- 
            - proj[, projectionAxes[1]]
    }
    if (flipYAxis) {
        proj[, projectionAxes[2]] <- 
            - proj[, projectionAxes[2]]
    }
    
    DF <- pData
    
    DF$x <- proj[, projectionAxes[1]]
    DF$y <- proj[, projectionAxes[2]]
    DF$stress <- mdsObj$spp
    
    xlabel <- paste0("Coord. ", projectionAxes[1])
    
    xlabel <- paste0(
        xlabel, 
        #" (marg. R2 : ", 
        #round(100*margRSq[projectionAxes[1]], 2), 
        " (% var. : ",
        round(100*explVar[projectionAxes[1]], 2),
        "%)")
    
    ylabel <- paste0("Coord. ", projectionAxes[2])
    
    ylabel <- paste0(
        ylabel, 
        #" (marg. R2 : ",
        #round(100*margRSq[projectionAxes[2]], 2), 
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
    
    mainAesMapping <- ggplot2::aes(
        x = .data[["x"]],
        y = .data[["y"]])
    
    maxAdditionalLabellingMapping <- 3
    if (!is.null(pDataForAdditionalLabelling) && !biplot) {
        if (!is.character(pDataForAdditionalLabelling)) {
            stop("pDataForAdditionalLabelling should be a character")
        }
        if (length(pDataForAdditionalLabelling) > maxAdditionalLabellingMapping) {
            stop("pDataForAdditionalLabelling length should be maximum 3")
        }
        if (length(pDataForAdditionalLabelling) >= 1) {
            mainAesMapping <- 
                c(mainAesMapping,
                  ggplot2::aes(text = .data[[pDataForAdditionalLabelling[1]]]))
        }
        if (length(pDataForAdditionalLabelling) >= 2) {
            mainAesMapping <- 
                c(mainAesMapping,
                  ggplot2::aes(text2 = .data[[pDataForAdditionalLabelling[2]]]))
        }
        if (length(pDataForAdditionalLabelling) >= 3) {
            mainAesMapping <- 
                c(mainAesMapping,
                  ggplot2::aes(text3 = .data[[pDataForAdditionalLabelling[3]]]))
        }
        
        # avoids error message: mapping should be created with `aes()`
        attr(mainAesMapping, "class") <- "uneval"
    }
    
    rangeMin <- min(proj)
    rangeMax <- max(proj)
    axesLimits <- c(rangeMin, rangeMax)
    #axesLimits <- c(min(-rangeMax, rangeMin), max(-rangeMin, rangeMax))
    
    ggplot2::scale_x_continuous(limits = axesLimits)
    
    p <- ggplot2::ggplot(
        data = DF,
        mapping = mainAesMapping) +
        ggplot2::labs(
            x = xlabel,
            y = ylabel,
            title = title,
            subtitle = subtitle) + 
        ggplot2::scale_x_continuous(limits = axesLimits) + 
        ggplot2::scale_y_continuous(limits = axesLimits)
    
    
    colourVar <- pDataForColour
    shapeVar <- pDataForShape
    
    geomPointMapping <- NULL
    
    if (!is.null(pDataForColour)) {
        if (!is.character(pDataForColour)) {
            stop("pDataForColour should be a character")
        }
        geomPointMapping <- c(
            geomPointMapping,
            ggplot2::aes(colour = .data[[colourVar]]))
    }
    if (!is.null(pDataForShape)) {
        if (!is.character(pDataForShape)) {
            stop("pDataForShape should be a character")
        }
        geomPointMapping <- c(
            geomPointMapping,
            ggplot2::aes(shape = .data[[shapeVar]]))
    }
    if (sizeReflectingStress) {
        geomPointMapping <- c(
            geomPointMapping,
            ggplot2::aes(size = .data[["stress"]]))
    }
    
    if (!is.null(geomPointMapping)) {
        # avoids error message: mapping should be created with `aes()`
        attr(geomPointMapping, "class") <- "uneval"
    }
    
    if (!is.null(pDataForLabel)) {
        if (!is.character(pDataForLabel)) {
            stop("pDataForLabel should be a character")
        }
        labelVar <- pDataForLabel
        if (repelPointsLabels) {
            p <- p + ggrepel::geom_text_repel(
                hjust=0.5, 
                vjust=1, 
                mapping = ggplot2::aes(label = .data[[labelVar]]),
                ...)
        } else {
            p <- p + ggplot2::geom_text(
                hjust=0.5, 
                vjust=1,
                mapping = ggplot2::aes(label = .data[[labelVar]]))
        }
    }
    
    p <- p + ggplot2::geom_point(mapping = geomPointMapping)
    
    # add biplot if specified
    if (biplot) {
        #browser()
        
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
            segmentXs <- rep(0., nExtVar)
            segmentYs <- rep(0., nExtVar)
            segmentNames <- rep("", nExtVar)
            visible <- rep(FALSE, nExtVar)
            for (j in seq_len(nExtVar)) {
                if (biplotType == "regression"){
                    segmentLength <- radius * mdsBiplot$R2vec[j]
                    segmentAngle <- atan(
                        mdsBiplot$coefficients[2,j] / 
                            mdsBiplot$coefficients[1,j])
                    if(mdsBiplot$coefficients[1,j] < 0){
                        segmentAngle <- segmentAngle + pi
                    }
                    segmentXs[j] <- segmentLength * cos(segmentAngle)
                    segmentYs[j] <- segmentLength * sin(segmentAngle)
                } else {
                    # biplotType == "correlation"
                    segmentXs[j] <- radius * mdsBiplot$pearsonCorr[1, j]
                    segmentYs[j] <- radius * mdsBiplot$pearsonCorr[2, j]
                    segmentLength <- sqrt(segmentXs[j]^2 + segmentYs[j]^2)
                }
                
                if (segmentLength >= lengthThreshold){
                    visible[j] <- TRUE
                    segmentNames[j] <- colnames(mdsBiplot$coefficients)[j]
                }
            }
            
            if (flipXAxis) {
                segmentXs <- - segmentXs   
            }
            if (flipYAxis) {
                segmentYs <- - segmentYs
            }
            
            segmentDF <- data.frame(
                segmentName = segmentNames,
                segmentXOrigin = rep(0., nExtVar),
                segmentYOrigin = rep(0., nExtVar),
                segmentX = segmentXs, 
                segmentY = segmentYs,
                visible = visible)
            
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
                    data = segmentDF)
                
            } else {
                # discarding possible warning message: 
                # 'Ignoring unknown aesthetics: text1'
                p <- p + ggplot2::geom_text(
                    mapping = newMapping,
                    data = segmentDF)
            }
            
            # add a dashed centered circle to obtain the RSq=1 benchmark
            p <- p + ggforce::geom_circle(
                mapping = ggplot2::aes(
                    x0 = 0.,
                    y0 = 0.,
                    r = radius),
                linetype = 'dotted')
        }
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
#' @param mdsObj a MDS object calculated by the SMACOF algorithm using
#' the computeMetricMDS() function
#' @param nDim number of dimensions to use when calculating   
#' Shepard's diagram and Rsquare.  
#' If `NULL`, it will be set equal to the number of projection dimensions  
#' as calculated in `mdsObj`
#' @param title title to give to the plot
#' @param pointSize plot size of points
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
#' # prepare data, build MDS object
#' example("computeMetricMDS") 
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
        nDim = NULL,
        title = "Multi Dimensional Scaling - Shepard's diagram",
        pointSize = 0.5,
        displayPseudoRSq = TRUE) {
    
    if (!inherits(mdsObj,"mdsRes")) {
        stop("mdsObj should be a 'mdsRes' object")
    }
    
    if (is.null(nDim)) {
        nDim <- ncol(mdsObj$proj)
    } else if (!is.numeric(nDim)) {
        stop("nDim should be numeric")
    } else if (nDim < 1) {
        stop("nDim should be >=1")
    } else if (ncol(mdsObj$proj) < nDim) {
        stop(
            "nDim too high compared to projection stored mdsObj (nDim = ",
            ncol(mdsObj$proj), ")")
    }
    
    RSq <- mdsObj$RSq[nDim]
    GoF <- mdsObj$GoF[nDim]
    
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
    
    projDist <- as.vector(dist(mdsObj$proj[,seq_len(nDim)]))
    HDDist <- as.vector(as.dist(mdsObj$pwDist))
    
    DF <- data.frame(
        HDDist = HDDist,
        projDist = projDist)
    
    p <- ggplot2::ggplot(
        data = DF,
        mapping = ggplot2::aes(x = .data[["HDDist"]],
                               y = .data[["projDist"]])) +
        
        ggplot2::labs(
            x = xlabel,
            y = ylabel,
            title = title,
            subtitle = subtitle) + 
        
        ggplot2::geom_point(colour = "blue", size = pointSize) +
        ggplot2::geom_abline(
            intercept = 0.,
            slope = 1.,
            linetype = 'dotted')
    p
}


#' @title SampleMDS biplot wrapping
#' @description `ggplotSampleMDSWrapBiplots` calls `ggplotSampleMDS` 
#' repeatly to generate biplots with different sets of external variables
#' and align them in a grid using the `patchwork` package, in a similar fashion 
#' as `ggplot2::facet_wrap()` does.
#' @param ncol passed to `patchwork::wrap_plots()`
#' @param nrow passed to `patchwork::wrap_plots()`
#' @param byrow passed to `patchwork::wrap_plots()`
#' @param mdsObj a MDS object calculated by the SMACOF algorithm using
#' the computeMetricMDS() function
#' @param pData a data.frame providing user input sample data. 
#' These can be design of experiment variables, phenotype data per sample,...
#' and will be used to highlight sample categories in the plot. 
#' @param projectionAxes which two axes should be plotted 
#' (should be a numeric vector of length 2)
#' @param extVariableList should be a named list of external variable matrices
#' Each element of the list should be a matrix with named columns 
#' corresponding to the variables. 
#' The number of rows should be the same as the number of samples. 
#' @param biplotType type of biplot used (see `ggplotSampleMDS()`)
#' @param pDataForColour if not NULL, which `pData` variable
#' will be used as colour aesthetic. Should be a character.
#' @param pDataForShape if not NULL, which `pData` variable
#' will be used as shape aesthetic. Should be a character.
#' @param pDataForLabel if not NULL, which `pData` variable 
#' will be used as point labels in the plot. Should be a character.
#' @param sizeReflectingStress if TRUE, size of points will appear 
#' proportional to stress by point, i.e. the bigger the sample point appears,
#' the less accurate its representation is 
#' (in terms of distances w.r.t. other points)
#' @param repelPointsLabels if TRUE, uses `ggrepel::geom_text_repel()` 
#' instead of `ggplot2::geom_text()`
#' (try to split the labels such that they do not overlap) for the points
#' @param repelArrowLabels if TRUE, uses `ggrepel::geom_text_repel()` 
#' instead of `ggplot2::geom_text()` for the arrows
#' @param displayArrowLabels if TRUE, displays arrows labels
#' @param arrowThreshold arrows will be made barely visible 
#' if their length is (in absolute value) less than this threshold.  
#' @param flipXAxis if TRUE, take the opposite of x values 
#' (provided as it might ease low dimensional projection comparisons)
#' @param flipYAxis if TRUE, take the opposite of y values 
#' (provided as it might ease low dimensional projection comparisons)
#' @param ... additional parameters passed to `ggrepel::geom_text_repel()` 
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
#' # prepare data, build MDS object
#' example("computeMetricMDS")
#' 
#' # plot mds projection on axes 1 and 2,
#' # use 'group' for colour, 'type' for shape, and no label 
#' 
#' p_12 <- ggplotSampleMDS(
#'     mdsObj = mdsObj,
#'     pData = flowCore::pData(fsAll),
#'     projectionAxes = c(1,2),
#'     pDataForColour = "grpId",
#'     pDataForLabel = NULL,
#'     pDataForShape = "type")
#' 
#' # try to associate axes with median or std deviation of each channel
#' # => use bi-plots
#' 
#' extVarList <- getChannelSummaryStats(
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
#'     pDataForLabel = NULL,
#'     pDataForShape = "type",
#'     seed = 0)
#' 
ggplotSampleMDSWrapBiplots <- function(
        mdsObj,
        extVariableList,
        ncol = NULL,
        nrow = NULL,
        byrow = NULL,
        pData = data.frame(
            sampleId = seq_len(nrow(mdsObj$proj))),
        projectionAxes = c(1,2),
        biplotType = c("correlation", "regression"),
        pDataForColour = NULL,
        pDataForShape = NULL,
        pDataForLabel = NULL,
        sizeReflectingStress = FALSE,
        repelPointsLabels = TRUE,
        displayArrowLabels = TRUE,
        arrowThreshold = 0.8,
        repelArrowLabels = FALSE,
        flipXAxis = FALSE,
        flipYAxis = FALSE,
        ...){
    
    #browser()
    
    if (!is.list(extVariableList)) {
        stop("[extVariableList] should be a non zero length list!")
    }
    
    if (is.null(names(extVariableList))) {
        stop("[extVariableList] should be a named list!")
    }
    
    nPlots <- length(extVariableList)
    
    pList <- list()
    
    for (i in seq_len(nPlots)) {
        p <- ggplotSampleMDS(
            mdsObj = mdsObj,
            pData = pData,
            title = paste0("biplot with ",
                           names(extVariableList)[i]),
            biplot = TRUE,
            biplotType = biplotType,
            extVariables = extVariableList[[i]],
            pDataForColour = pDataForColour,
            pDataForShape = pDataForShape,
            pDataForLabel = pDataForLabel,
            pDataForAdditionalLabelling = NULL,
            sizeReflectingStress = sizeReflectingStress,
            repelPointsLabels = repelPointsLabels,
            displayArrowLabels = displayArrowLabels,
            arrowThreshold = arrowThreshold,
            repelArrowLabels = repelArrowLabels,
            flipXAxis = flipXAxis,
            flipYAxis = flipYAxis,
            ...) 
        pList[[i]] <- p + ggplot2::labs(subtitle = NULL)
    }
        
    p <- patchwork::wrap_plots(
        pList, 
        ncol = ncol,
        nrow = nrow,
        byrow = byrow,
        guides = 'collect')
    p
}


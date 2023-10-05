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
#' @description `ggplotSamplesMDS` uses ggplot2 
#' to provide plots of Metric MDS results.   
#' If both `projectionAxes` are in c(1,2), 
#' a global RSquare for these 2 axes is provided. For additional dimensions,
#' a marginal RSquare per dimension is provided.
#' @rdname CytoMDS
#' @param mdsObj a MDS object calculated by the SMACOF algorithm using
#' the computeMetricMDS() function
#' @param pData a data.frame providing user input sample data. 
#' These can be design of experiment variables, phenotype data per sample,...
#' and will be used to highlight sample categories in the plot. 
#' @param projectionAxes which two axes should be plotted 
#' (should be a numeric vector of length 2)
#' @param biplot if TRUE, adds projection of external variables
#' @param extVariables are used to generate a biplot
#' these are the external variables to regress with obtained configuration
#' according to the two projection axes
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
#' @param repelArrowsLabels if TRUE, uses `ggrepel::geom_text_repel()` 
#' instead of `ggplot2::geom_text()` for the arrows
#' @param flipXAxis if TRUE, take the opposite of x values 
#' (provided as it might ease low dimensional projection comparisons)
#' @param flipYAxis if TRUE, take the opposite of y values 
#' (provided as it might ease low dimensional projection comparisons)
#' @param ... additional parameters passed to `ggrepel::geom_text_repel()` 
#' (if used)
#' @importFrom stats as.dist dist lm 
#' @importFrom rlang .data
#' 
#' @export
#' 
ggplotSamplesMDS <- function(
        mdsObj,
        pData = data.frame(
            sampleId = seq_len(nrow(mdsObj$proj))),
        projectionAxes = c(1,2),
        biplot = FALSE,
        extVariables = NULL,
        pDataForColour = NULL,
        pDataForShape = NULL,
        pDataForLabel = "name",
        pDataForAdditionalLabelling = NULL,
        sizeReflectingStress = FALSE,
        title = "Multi Dimensional Scaling",
        repelPointsLabels = TRUE,
        repelArrowsLabels = FALSE,
        flipXAxis = FALSE,
        flipYAxis = FALSE,
        ...){
        
    #browser()
    
    if (!inherits(mdsObj, "mdsRes")) {
        stop("mdsObj should be a 'mdsRes' object")
    }
    
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

    subtitle <- paste0("(Pseudo R2 = ", round(RSq, 4),
                       #"; Goodness of Fit = ", round(GoF, 4),
                       "; nDim = ", nDim, ")")
    
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
        
        mdsBiplot <- computeMetricMDSBiplot(
            mdsObj,
            projectionAxes = projectionAxes,
            extVariables = extVariables)
            
        radius <- min(-axesLimits[1], axesLimits[2])
        nExtVar <- ncol(mdsBiplot$coefficients)
        if (nExtVar > 0) {
            segmentXs <- rep(0., nExtVar)
            segmentYs <- rep(0., nExtVar)
            segmentNames <- rep("", nExtVar)
            for (j in seq_len(nExtVar)) {
                segmentLength <- radius * mdsBiplot$R2vec[j]
                segmentAngle <- atan(
                    mdsBiplot$coefficients[2,j] / 
                        mdsBiplot$coefficients[1,j])
                if(mdsBiplot$coefficients[1,j] < 0){
                    segmentAngle <- segmentAngle + pi
                }
                segmentXs[j] <- segmentLength * cos(segmentAngle)
                segmentYs[j] <- segmentLength * sin(segmentAngle)
                segmentNames[j] <- colnames(mdsBiplot$coefficients)[j]
            }
            
            segmentDF <- data.frame(
                segmentName = segmentNames,
                segmentXOrigin = rep(0., nExtVar),
                segmentYOrigin = rep(0., nExtVar),
                segmentX = segmentXs, 
                segmentY = segmentYs)
            
            p <- p + ggplot2::geom_segment(
                mapping = ggplot2::aes(
                    x = .data[["segmentXOrigin"]],
                    y = .data[["segmentYOrigin"]],
                    xend = .data[["segmentX"]],
                    yend = .data[["segmentY"]]),
                data = segmentDF,
                arrow = ggplot2::arrow(
                    length = ggplot2::unit(0.1, "inches"),
                    type = "closed"),
                linetype = 'dotted')
            
            
            if (repelArrowsLabels) {
                p <- p + ggrepel::geom_text_repel(
                    mapping = ggplot2::aes(
                        x = (.data[["segmentXOrigin"]] + 
                                 .data[["segmentX"]]) / 2,
                        y = (.data[["segmentYOrigin"]] + 
                                 .data[["segmentY"]]) / 2, 
                        label = .data[["segmentName"]]),
                    data = segmentDF)
                
            } else {
                p <- p + ggplot2::geom_text(
                    mapping = ggplot2::aes(
                        x = (.data[["segmentXOrigin"]] + 
                                 .data[["segmentX"]]) / 2,
                        y = (.data[["segmentYOrigin"]] + 
                                 .data[["segmentY"]]) / 2, 
                        label = .data[["segmentName"]]),
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
#' @description `ggplotSamplesMDSShepard` uses ggplot2 
#' to provide plot of Metric MDS results.  
#' Shepard diagram provides a scatter plot of :
#' - on the x axis, the high dimensional pairwise distances 
#' between each sample pairs  
#' - on the y axis, the corresponding pairwise distances in the obtained 
#' low dimensional projection
#' @rdname CytoMDS
#' @param mdsObj a MDS object calculated by the SMACOF algorithm using
#' the computeMetricMDS() function
#' @param nDim number of dimensions to use when calculating   
#' Shepard's diagram and Rsquare.  
#' If `NULL`, it will be set equal to the number of projection dimensions  
#' as calculated in `mdsObj`
#' @param title title to give to the plot
#' @param pointSize plot size of points
#'
#' @importFrom rlang .data
#' @export
#' 
ggplotSamplesMDSShepard <- function(
        mdsObj,
        nDim = NULL,
        title = "Multi Dimensional Scaling - Shepard's diagram",
        pointSize = 0.5) {
    
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
    subtitle <- paste0("(Pseudo R2 = ", round(RSq, 4),
                       #"; Goodness of Fit = ", round(GoF, 4),
                       "; nDim = ", nDim, ")")
    
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



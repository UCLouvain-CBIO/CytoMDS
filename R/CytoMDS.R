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


#' @title CytoMDS package
#' 
#' @name CytoMDS
#'
#' @rdname CytoMDS
#' 
#' @seealso [ggplotSamplesMDS], [ggplotSamplesMDSShepard], [computeMetricMDS]
#'
#' @description
#'
#' `CytoMDS` implements a low dimensional visualization of a set of cytometry 
#' samples, in order to visually assess the 'distances' between them.
#' This, in turn, can greatly help the user to identify quality issues 
#' like batch effects or outlier samples, and/or check the presence of 
#' potential sample clusters that might align with the experimental design.  
#'   
#' The CytoMDS algorithm combines, on the one hand, 
#' the concept of Earth Mover's Distance (EMD), a.k.a. Wasserstein metric 
#' and, on the other hand, the Multi Dimensional Scaling (MDS) algorithm 
#' for the low dimensional projection.  
#' 
#' Also, the package provides some diagnostic tools for 
#' both checking the quality of the MDS projection, 
#' as well as tools to help with the interpretation of 
#' the axes of the projection.
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
#' pwDist <- getPairWiseEMDDist(fsAll, 
#'                              channels = c("FSC-A", "SSC-A"),
#'                              verbose = FALSE)
#' 
#' # compute Metric MDS object
#' 
#' mdsObj <- computeMetricMDS(pwDist, nDim = 4, seed = 0)
#' 
#' # plot mds projection on axes 1 and 2,
#' # use 'group' for colour, 'type' for shape, and no label 
#' 
#' p_12 <- ggplotSamplesMDS(
#'     mdsObj = mdsObj,
#'     pData = flowCore::pData(fsAll),
#'     projectionAxes = c(1,2),
#'     pDataForColour = "group",
#'     pDataForLabel = NULL,
#'     pDataForShape = "type")
#' 
#' # plot mds projection on axes 3 and 4,
#' # use 'group' for colour, and 'name' as point label
#' 
#' p_34 <- ggplotSamplesMDS(
#'     mdsObj = mdsObj,
#'     pData = flowCore::pData(fsAll),
#'     projectionAxes = c(3,4),
#'     pDataForColour = "group",
#'     pDataForLabel = "name")
#' 
#' # plot mds projection on axes 1 and 2,
#' # use 'group' for colour, 'type' for shape, and 'name' as point label
#' # have sample point size reflecting 'stress'
#' # i.e. quality of projection w.r.t. distances to other points
#' 
#' p12_Stress <- ggplotSamplesMDS(
#'     mdsObj = mdsObj,
#'     pData = flowCore::pData(fsAll),
#'     projectionAxes = c(1,2),
#'     pDataForColour = "group",
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
#' bp_12 <- ggplotSamplesMDS(
#'     mdsObj = mdsObj,
#'     pData = flowCore::pData(fsAll),
#'     projectionAxes = c(1,2),
#'     biplot = TRUE,
#'     extVariables = extVars,
#'     pDataForColour = "group",
#'     pDataForLabel = NULL,
#'     pDataForShape = "type",
#'     seed = 0)
#' 
#' bp_34 <- ggplotSamplesMDS(
#'     mdsObj = mdsObj,
#'     pData = flowCore::pData(fsAll),
#'     projectionAxes = c(3,4),
#'     biplot = TRUE,
#'     extVariables = extVars,
#'     pDataForColour = "group",
#'     pDataForLabel = "name",
#'     seed = 0)
#' 
#' # Shepard diagrams 
#' 
#' p2D <- ggplotSamplesMDSShepard(
#'     mdsObj,
#'     nDim = 2,
#'     pointSize = 1,
#'     title = "Shepard with 2 dimensions")
#' 
#' p3D <- ggplotSamplesMDSShepard(
#'     mdsObj,
#'     nDim = 3,
#'     title = "Shepard with 3 dimensions") 
#'     #' 
#' pDefD <- ggplotSamplesMDSShepard(
#'     mdsObj,
#'     title = "Shepard with default nb of dimensions") 
#'
NULL
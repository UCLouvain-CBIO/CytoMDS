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

# loading test dataset from CytoPipeline
library(CytoPipeline)

data(OMIP021Samples)

transList <- estimateScaleTransforms(
    ff = OMIP021Samples[[1]],
    fluoMethod = "estimateLogicle",
    scatterMethod = "linearQuantile",
    scatterRefMarker = "BV785 - CD3")

OMIP021Trans <- CytoPipeline::applyScaleTransforms(
    OMIP021Samples, 
    transList)

ffList <- flowCore::flowSet_to_list(OMIP021Trans)

for(i in 3:5){
    ffList[[i]] <- 
        aggregateAndSample(
            OMIP021Trans,
            seed = 10*i,
            nTotalEvents = 5000)[,1:22]
}

fsNames <- c("Donor1", "Donor2", paste0("Agg",1:3))
names(ffList) <- fsNames

fsAll <- as(ffList,"flowSet")

flowCore::pData(fsAll)$type <- factor(c("real", "real", rep("synthetic", 3)))
flowCore::pData(fsAll)$grpId <- factor(c("D1", "D2", rep("Agg", 3)))

pwDist <- pairwiseEMDDist(fsAll, 
                             channels = c("FSC-A", "SSC-A"),
                             verbose = FALSE)

test_that("ggplotSampleMDS works", {
    
    mdsObj <- computeMetricMDS(pwDist, nDim = 4, seed = 0)
    
    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(1,2),
                         pDataForColour = "grpId",
                         pDataForLabel = NULL,
                         pDataForShape = "type")
    
    vdiffr::expect_doppelganger("ggplotSampleMDS with axes 1 and 2",
                                fig = p)
    
    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(3,4),
                         pDataForColour = "grpId",
                         pDataForLabel = "name",
                         seed = 0)
    
    vdiffr::expect_doppelganger(
        "ggplotSampleMDS with axes 3 and 4",
        fig = p)
    
    extVars <- getChannelSummaryStats(
        fsAll,
        channels = c("FSC-A", "SSC-A"),
        statFUNs = stats::median,
        verbose = FALSE)
    
    
    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(1,2),
                         biplot = TRUE,
                         extVariables = extVars,
                         pDataForColour = "grpId",
                         pDataForLabel = NULL,
                         pDataForShape = "type",
                         seed = 0)
    
    vdiffr::expect_doppelganger(
        "ggplotSampleMDS with axes 1 and 2 and extVars",
        fig = p)
    
    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(1,2),
                         biplot = TRUE,
                         biplotType = "regression",
                         extVariables = extVars,
                         pDataForColour = "grpId",
                         pDataForLabel = NULL,
                         pDataForShape = "type",
                         seed = 0)
    
    vdiffr::expect_doppelganger(
        "ggplotSampleMDS axes 1 2 biplot regression",
        fig = p)    
    
    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(3,4),
                         biplot = TRUE,
                         extVariables = extVars,
                         pDataForColour = "grpId",
                         pDataForLabel = "name",
                         seed = 0)
    
    vdiffr::expect_doppelganger(
        "ggplotSampleMDS with axes 3 and 4 and extVars",
        fig = p)
    
    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(3,4),
                         biplot = TRUE,
                         biplotType = "regression",
                         extVariables = extVars,
                         pDataForColour = "grpId",
                         pDataForLabel = NULL,
                         pDataForShape = "type",
                         seed = 0)
    
    vdiffr::expect_doppelganger(
        "ggplotSampleMDS axes 3 4 biplot regression",
        fig = p)    
    
    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(3,4),
                         biplot = TRUE,
                         extVariables = extVars,
                         pDataForColour = "grpId",
                         pDataForLabel = "name",
                         arrowThreshold = 0.,
                         seed = 0)
    
    vdiffr::expect_doppelganger(
        "ggplotSampleMDS arrowThreshold",
        fig = p)
    
    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(1,2),
                         biplot = TRUE,
                         extVariables = extVars,
                         pDataForColour = "grpId",
                         pDataForLabel = NULL,
                         pDataForShape = "type",
                         displayArrowLabels = FALSE,
                         seed = 0)
    
    vdiffr::expect_doppelganger(
        "ggplotSampleMDS no arrow label",
        fig = p)
    
    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(3,4),
                         biplot = TRUE,
                         extVariables = extVars,
                         pDataForColour = "grpId",
                         pDataForLabel = NULL,
                         seed = 0)
    
    vdiffr::expect_doppelganger(
        "ggplotSampleMDS no point labels",
        fig = p)
    
    # use 2 dimensions to make stress per point meaningful
    mdsObj <- computeMetricMDS(pwDist, nDim = 2, seed = 0)
    
    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(1,2),
                         pDataForColour = "grpId",
                         pDataForLabel = "name",
                         pDataForShape = "type",
                         sizeReflectingStress = TRUE,
                         seed = 0)
    
    vdiffr::expect_doppelganger("ggplotSampleMDS with sizeReflectingStress",
                                fig = p)
    
    # use pData for additional labelling 
    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(1,2),
                         pDataForAdditionalLabelling = c("grpId", "type"),
                         repelPointsLabels = FALSE)
    
    expect_equal(p$labels$text, "grpId")
    expect_equal(p$labels$text2, "type")
    
    # test that pData for additional labelling has been well ignored
    # with biplot activated
    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         biplot = TRUE,
                         extVariables = extVars,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(1,2),
                         pDataForAdditionalLabelling = c("grpId", "type"),
                         repelPointsLabels = FALSE)
    
    expect_null(p$labels$text)
    expect_null(p$labels$text2)
    
    # test flipXAxis and flipYAxis
    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(1,2),
                         pDataForColour = "grpId",
                         pDataForLabel = "name",
                         pDataForShape = "type",
                         seed = 0,
                         flipXAxis = TRUE)
    
    vdiffr::expect_doppelganger("ggplotSampleMDS with flipX",
                                fig = p)
    
    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(1,2),
                         pDataForColour = "grpId",
                         pDataForLabel = "name",
                         pDataForShape = "type",
                         seed = 0,
                         flipXAxis = TRUE,
                         flipYAxis = TRUE)
    
    vdiffr::expect_doppelganger("ggplotSampleMDS with flipX-Y",
                                fig = p)
    
    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         biplot = TRUE,
                         extVariables = extVars,
                         projectionAxes = c(1,2),
                         pDataForAdditionalLabelling = c("grpId", "type"),
                         repelPointsLabels = FALSE,
                         flipXAxis = TRUE,
                         flipYAxis = TRUE)
    
    vdiffr::expect_doppelganger("ggplotSampMDS with bipl-flpX-Y",
                                fig = p)
    
})


test_that("ggplotSampleMDSWrapBiplots works", {
    mdsObj <- computeMetricMDS(pwDist, nDim = 4, seed = 0)
    
    # try to associate axes with median or std deviation of each channel
    # => use bi-plots

    extVarList <- getChannelSummaryStats(
        fsAll,
        channels = c("FSC-A", "SSC-A"),
        statFUNs = c("median" = stats::median, 
                     "std.dev" = stats::sd))

    bpFull <- ggplotSampleMDSWrapBiplots(
        mdsObj = mdsObj,
        extVariableList = extVarList,
        pData = flowCore::pData(fsAll),
        projectionAxes = c(1,2),
        pDataForColour = "grpId",
        pDataForLabel = NULL,
        pDataForShape = "type",
        seed = 0)
    
    vdiffr::expect_doppelganger(
        "ggplotSampleMDSWrapBiplots default rows and cols",
        fig = bpFull)
    
    bpFull2 <- ggplotSampleMDSWrapBiplots(
        mdsObj = mdsObj,
        extVariableList = extVarList,
        ncol = 1,
        pData = flowCore::pData(fsAll),
        projectionAxes = c(1,2),
        pDataForColour = "grpId",
        pDataForLabel = NULL,
        pDataForShape = "type",
        seed = 0)
    
    vdiffr::expect_doppelganger(
        "ggplotSampleMDSWrapBiplots with 1 col",
        fig = bpFull2)
})

test_that("ggplotSampleMDSShepard works", {
    
    mdsObj <- computeMetricMDS(pwDist, nDim = 4, seed = 0)
    
    p <- ggplotSampleMDSShepard(mdsObj,
                                nDim = 2,
                                pointSize = 1,
                                title = "Shepard with 2 dimensions")
    
    vdiffr::expect_doppelganger(
        "ggplotSampleMDSShepard with 2 dimensions",
        fig = p)
    
    p <- ggplotSampleMDSShepard(mdsObj,
                                nDim = 3,
                                title = "Shepard with 3 dimensions") 
    
    vdiffr::expect_doppelganger(
        "ggplotSampleMDSShepard with 3 dimensions",
        fig = p)
    
    p <- ggplotSampleMDSShepard(mdsObj,
                                title = "Shepard with default nb of dimensions") 
    
    vdiffr::expect_doppelganger(
        "ggplotSampleMDSShepard with default dim nb",
        fig = p)
    
})

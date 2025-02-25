# CytoMDS - Copyright (C) <2023-2024>
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

# work-around for temporary fail on GHA of some graphical tests
activate_ggplotWrapPlots <- FALSE

data(OMIP021Samples)

transList <- estimateScaleTransforms(
    ff = OMIP021Samples[[1]],
    fluoMethod = "estimateLogicle",
    scatterMethod = "linearQuantile",
    scatterRefMarker = "BV785 - CD3")

OMIP021Trans <- CytoPipeline::applyScaleTransforms(
    OMIP021Samples, 
    transList)

ffList <- c(
    flowCore::flowSet_to_list(OMIP021Trans),
    lapply(3:5,
           FUN = function(i) {
               aggregateAndSample(
                   OMIP021Trans,
                   seed = 10*i,
                   nTotalEvents = 5000)[,1:22]
           }))

fsNames <- c("Donor1", "Donor2", paste0("Agg",1:3))
names(ffList) <- fsNames

fsAll <- as(ffList,"flowSet")

fsAll <- as(ffList,"flowSet")

flowCore::pData(fsAll)$type <- factor(c("real", "real", rep("synthetic", 3)))
flowCore::pData(fsAll)$grpId <- factor(c("D1", "D2", rep("Agg", 3)))
flowCore::pData(fsAll)$lbl <- paste0("S", 1:5)

pwDist <- pairwiseEMDDist(fsAll, 
                          channels = c("FSC-A", "SSC-A"),
                          verbose = FALSE)

test_that("ggplotMarginalDensities works", {
    
    p <- ggplotMarginalDensities(
        OMIP021Samples,
        pDataForGroup = "Donor",
        pDataForColour = "Donor",
        transList = transList
    )
    vdiffr::expect_doppelganger(
        "ggplotMarginalDensities no channels with transList",
        fig = p)
    
    selChannels <- c("FSC-A", "SSC-A", "670/30Violet-A", "525/50Violet-A")
    p <- ggplotMarginalDensities(
        OMIP021Samples,
        channels = selChannels,
        pDataForGroup = "Donor",
        pDataForColour = "Donor",
        transList = transList
    )
    vdiffr::expect_doppelganger(
        "ggplotMarginalDensities with channels with transList",
        fig = p)
    
    selChannels <- c("FSC-A", "SSC-A", "BV785 - CD3", "APCCy7 - CD4")
    p <- ggplotMarginalDensities(
        OMIP021Samples,
        channels = selChannels,
        pDataForGroup = "Donor",
        pDataForColour = "Donor",
        transList = transList
    )
    vdiffr::expect_doppelganger(
        "ggplotMarginalDensities with markers with transList",
        fig = p)
    
    p <- ggplotMarginalDensities(
        fsAll
    )
    vdiffr::expect_doppelganger(
        "ggplotMarginalDensities nothing",
        fig = p)
    
    
    p <- ggplotMarginalDensities(
        fsAll,
        pDataForGroup = "lbl"
    )
    vdiffr::expect_doppelganger(
        "ggplotMarginalDensities with groupBy",
        fig = p)
    
    p <- ggplotMarginalDensities(
        fsAll,
        pDataForGroup = "lbl",
        pDataForColour = "grpId"
    )
    vdiffr::expect_doppelganger(
        "ggplotMarginalDensities with groupBy and colourBy",
        fig = p)
    
    p <- ggplotMarginalDensities(
        fsAll,
        pDataForGroup = "lbl",
        pDataForColour = "lbl"
    )
    vdiffr::expect_doppelganger(
        "ggplotMarginalDensities with groupBy and same colourBy",
        fig = p)
    
    # sample subset
    
    p <- ggplotMarginalDensities(
        fsAll,
        sampleSubset = pData(phenoData(fsAll))[, "type"] == "synthetic",
        pDataForGroup = "lbl",
        pDataForColour = "lbl"
    )
    
    vdiffr::expect_doppelganger(
        "ggplotMarginalDensities with sample subset",
        fig = p)
    
    # subsampling
    
    p <- ggplotMarginalDensities(
        fsAll,
        nEventInSubsample = 100,
        seed = 0,
        pDataForGroup = "lbl",
        pDataForColour = "lbl"
    )
    
    vdiffr::expect_doppelganger(
        "ggplotMarginalDensities with subsampling",
        fig = p)
    
    # wrong input class
    expect_error(ggplotMarginalDensities("hello"),
                 regexp = "should be a flowSet or a flowFrame")
    
    # wrong pDataForColour
    expect_error(ggplotMarginalDensities(
        fsAll,
        pDataForColour = "labbbbbbel"),
        regexp = "'pDataForColour' should be in phenoData columns"
    )
    
    # wrong pDataForGroup
    expect_error(ggplotMarginalDensities(
        fsAll,
        pDataForGroup = "grrrrroup"),
        regexp = "'pDataForGroup' should be in phenoData columns"
    )
    
    # with flowFrame
    p <- ggplotMarginalDensities(
        OMIP021Trans[[1]]
    )
    
    vdiffr::expect_doppelganger(
        "ggplotMarginalDensities with flowFrame",
        fig = p)
})

test_that("ggplotDistFeatureImportance works", {
    p <- ggplotDistFeatureImportance(pwDist)
    vdiffr::expect_doppelganger("ggplotDistFeatureImportance",
                                fig = p)
})

test_that("ggplotSampleMDS works", {

    mdsObj <- computeMetricMDS(pwDist, nDim = 4, seed = 0)

    set.seed(0) # to get same results with ggrepel()

    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(1,2),
                         pDataForColour = "grpId",
                         pDataForShape = "type")

    vdiffr::expect_doppelganger("ggplotSampleMDS with axes 1 and 2",
                                fig = p)

    # no labels

    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(1,2),
                         pDataForColour = "grpId",
                         pDataForShape = "type",
                         displayPointLabels = FALSE)

    vdiffr::expect_doppelganger("ggplotSampleMDS with axes 1 and 2 - no labels",
                                fig = p)

    # explicit labels

    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(1,2),
                         pDataForColour = "grpId",
                         pDataForShape = "type",
                         pDataForLabel = "lbl")

    vdiffr::expect_doppelganger(
        "ggplotSampleMDS with axes 1 and 2 - explicit labels",
        fig = p)

    # testing with subset

    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         sampleSubset = (flowCore::pData(fsAll)$type == "real"),
                         projectionAxes = c(1,2),
                         pDataForColour = "grpId",
                         pDataForShape = "type")

    vdiffr::expect_doppelganger("ggplotSampleMDS with axes 1 and 2 - real",
                                fig = p)

    expect_error(ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         sampleSubset = (grpId == "Agg"),
                         projectionAxes = c(1,2),
                         pDataForColour = "grpId",
                         pDataForShape = "type"),
                 regexp = "object 'grpId' not found")


    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(3,4),
                         pDataForColour = "grpId",
                         pDataForLabel = "name",
                         seed = 0)

    vdiffr::expect_doppelganger(
        "ggplotSampleMDS with axes 3 and 4",
        fig = p)

    extVars <- channelSummaryStats(
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
                         pDataForShape = "type",
                         seed = 0)

    vdiffr::expect_doppelganger(
        "ggplotSampleMDS axes 1 2 biplot regression",
        fig = p)
    
    # test arrow label size
    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(1,2),
                         biplot = TRUE,
                         biplotType = "regression",
                         extVariables = extVars,
                         pDataForColour = "grpId",
                         pDataForShape = "type",
                         arrowLabelSize = 6,
                         seed = 0)
    
    vdiffr::expect_doppelganger(
        "ggplotSampleMDS axes 1 2 biplot arrow label size",
        fig = p)

    extVarNAs <-  extVars
    extVarNAs[3,1] <- NA

    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(1,2),
                         biplot = TRUE,
                         extVariables = extVarNAs,
                         pDataForColour = "grpId",
                         pDataForShape = "type",
                         seed = 0)

    vdiffr::expect_doppelganger(
        "ggplotSampleMDS with axes 1 and 2 and extVars nas",
        fig = p)

    extVarInvalid <-  extVars
    extVarInvalid[,2] <- rep(5, extVars[1,2])

    expect_warning(pI <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(1,2),
                         biplot = TRUE,
                         extVariables = extVarInvalid,
                         pDataForColour = "grpId",
                         pDataForShape = "type",
                         seed = 0),
                   regexp = "external variable SSC-A is constant => discarded")

    vdiffr::expect_doppelganger(
        "ggplotSampleMDS with axes 1 and 2 and extVars invalid",
        fig = pI)

    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(1,2),
                         biplot = TRUE,
                         biplotType = "regression",
                         extVariables = extVarNAs,
                         pDataForColour = "grpId",
                         pDataForShape = "type",
                         seed = 0)

    vdiffr::expect_doppelganger(
        "ggplotSampleMDS axes 1 2 biplot regression nas",
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

    # testing with subset

    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         sampleSubset = (flowCore::pData(fsAll)$grpId == "Agg"),
                         projectionAxes = c(1,2),
                         biplot = TRUE,
                         extVariables = extVars,
                         pDataForColour = "grpId",
                         pDataForLabel = "name",
                         arrowThreshold = 0.,
                         seed = 0)

    vdiffr::expect_doppelganger(
        "ggplotSampleMDS arrowThreshold subset",
        fig = p)

    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(1,2),
                         biplot = TRUE,
                         extVariables = extVars,
                         pDataForColour = "grpId",
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
                         pointSizeReflectingStress = TRUE,
                         seed = 0)

    vdiffr::expect_doppelganger("ggplotSampleMDS with sizeReflectingStress",
                                fig = p)
    
    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(1,2),
                         pDataForColour = "grpId",
                         pDataForLabel = "name",
                         pDataForShape = "type",
                         pointSize = 2,
                         seed = 0)

    vdiffr::expect_doppelganger("ggplotSampleMDS with pointSize",
                                fig = p)

    # use pData for additional labelling
    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(1,2),
                         pDataForAdditionalLabelling = c("grpId", "type"),
                         repelPointLabels = FALSE)

    expect_equal(p$labels$text2, "grpId")
    expect_equal(p$labels$text3, "type")

    # test that pData for additional labelling has been well ignored
    # with biplot activated
    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         biplot = TRUE,
                         extVariables = extVars,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(1,2),
                         pDataForAdditionalLabelling = c("grpId", "type"),
                         repelPointLabels = FALSE)

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
                         repelPointLabels = FALSE,
                         flipXAxis = TRUE,
                         flipYAxis = TRUE)

    vdiffr::expect_doppelganger("ggplotSampMDS with bipl-flpX-Y",
                                fig = p)

    # test minimal call without pData argument
    p <- ggplotSampleMDS(mdsObj = mdsObj)

    vdiffr::expect_doppelganger("ggplotSampleMDS minimal call",
                                fig = p)
    
    # test point label size 
    p <- ggplotSampleMDS(mdsObj = mdsObj,
                         pData = flowCore::pData(fsAll),
                         projectionAxes = c(1,2),
                         pDataForColour = "grpId",
                         pDataForLabel = "name",
                         pDataForShape = "type",
                         seed = 0,
                         pointLabelSize = 6)
    
    vdiffr::expect_doppelganger("ggplotSampleMDS with pointLabelSize",
                                fig = p)

})


test_that("ggplotSampleMDSWrapBiplots works", {
    mdsObj <- computeMetricMDS(pwDist, nDim = 4, seed = 0)

    # try to associate axes with median or std deviation of each channel
    # => use bi-plots

    extVarList <- channelSummaryStats(
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
        pDataForShape = "type",
        seed = 0)

    if (activate_ggplotWrapPlots) {
        vdiffr::expect_doppelganger(
            "ggplotSampleMDSWrapBiplots default rows and cols",
            fig = bpFull)
    }
    
    bpFull <- ggplotSampleMDSWrapBiplots(
        mdsObj = mdsObj,
        extVariableList = extVarList,
        displayLegend = FALSE,
        pData = flowCore::pData(fsAll),
        projectionAxes = c(1,2),
        pDataForColour = "grpId",
        pDataForShape = "type",
        seed = 0)

    if (activate_ggplotWrapPlots) {
        vdiffr::expect_doppelganger(
            "ggplotSampleMDSWrapBiplots no legend",
            fig = bpFull)
    }
    
    # with subset
    bpFull <- ggplotSampleMDSWrapBiplots(
        mdsObj = mdsObj,
        extVariableList = extVarList,
        pData = flowCore::pData(fsAll),
        sampleSubset = (flowCore::pData(fsAll)$type == "synthetic"),
        projectionAxes = c(1,2),
        pDataForColour = "grpId",
        pDataForShape = "type",
        seed = 0)

    if (activate_ggplotWrapPlots) {
        vdiffr::expect_doppelganger(
        "ggplotSampleMDSWrapBiplots default rows and cols - subset",
        fig = bpFull)
    }

    bpFull2 <- ggplotSampleMDSWrapBiplots(
        mdsObj = mdsObj,
        extVariableList = extVarList,
        ncol = 1,
        pData = flowCore::pData(fsAll),
        projectionAxes = c(1,2),
        pDataForColour = "grpId",
        pDataForShape = "type",
        seed = 0)

    if (activate_ggplotWrapPlots) {
        vdiffr::expect_doppelganger(
            "ggplotSampleMDSWrapBiplots with 1 col",
            fig = bpFull2)
    }
    
    dum <- 0
    expect_equal(dum, 0)
})

test_that("ggplotSampleMDSShepard works", {

    set.seed(0) # to get same results with ggrepel()

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
    
    p <- ggplotSampleMDSShepard(mdsObj,
                                pointSize = 1,
                                lineWidth = 1)
    
    vdiffr::expect_doppelganger(
        "ggplotSampleMDSShepard with explicit graphical params",
        fig = p)
    

})

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

test_that("ggplotSamplesMDS works", {
    transList <- estimateScaleTransforms(
        ff = OMIP021Samples[[1]],
        fluoMethod = "estimateLogicle",
        scatterMethod = "linearQuantile",
        scatterRefMarker = "BV785 - CD3")
    
    ffList <- flowCore::flowSet_to_list(OMIP021Samples)
    
    for(i in 3:5){
        ffList[[i]] <- 
            aggregateAndSample(
                OMIP021Samples,
                seed = 10*i,
                nTotalEvents = 5000)[,1:22]
    }
    
    fsNames <- c("Donor1", "Donor2", paste0("Agg",1:3))
    names(ffList) <- fsNames
    
    fsAll <- as(ffList,"flowSet")
    
    flowCore::pData(fsAll)$group <- factor(c("D1", "D2", rep("Agg", 3)))
    flowCore::pData(fsAll)$type <- factor(c("real", "real", rep("synthetic", 3)))
    
    pwDist <- getPairWiseEMDDist(fsAll, 
                                 channels = c("FSC-A", "SSC-A"),
                                 transList = transList,
                                 verbose = FALSE)
    
    
    
    
    mdsObj <- computeMetricMDS(pwDist, nDim = 4, seed = 0)
    
    expect_equal(mdsObj$stress, 4.322212e-16)
    tgtspp <- c(
        13.5257732, 47.5051546, 12.6391753, 13.0309278, 13.2989691)
    names(tgtspp) <- 1:5
    expect_equal(mdsObj$spp, tgtspp)
    
    
    p <- ggplotSamplesMDS(mdsObj = mdsObj,
                          pData = flowCore::pData(fsAll),
                          projectionAxes = c(1,2),
                          pDataForColour = "group",
                          pDataForLabel = NULL,
                          pDataForShape = "type")
    
    vdiffr::expect_doppelganger("ggplotSamplesMDS with axes 1 and 2",
                                fig = p)
    
    p <- ggplotSamplesMDS(mdsObj = mdsObj,
                          pData = flowCore::pData(fsAll),
                          projectionAxes = c(3,4),
                          pDataForColour = "group",
                          pDataForLabel = "name",
                          seed = 0)
    
    vdiffr::expect_doppelganger(
        "ggplotSamplesMDS with axes 3 and 4",
        fig = p)
    
    p <- ggplotSamplesMDS(mdsObj = mdsObj,
                          pData = flowCore::pData(fsAll),
                          projectionAxes = c(1,2),
                          pDataForColour = "group",
                          pDataForLabel = "name",
                          pDataForShape = "type",
                          sizeReflectingStress = TRUE,
                          seed = 0)
    
    vdiffr::expect_doppelganger("ggplotSamplesMDS with sizeReflectingStress",
                                fig = p)
    
    extVars <- getChannelsSummaryStat(
        fsAll,
        channels = c("FSC-A", "SSC-A"),
        transList = transList,
        statFUN = stats::median,
        verbose = FALSE)
    
    # mdsObj <- computeMetricMDS(pwDist, nDim = 4, seed = 0,
    #                            extVariables = extVars)
    
    
    p <- ggplotSamplesMDS(mdsObj = mdsObj,
                          pData = flowCore::pData(fsAll),
                          projectionAxes = c(1,2),
                          biplot = TRUE,
                          extVariables = extVars,
                          pDataForColour = "group",
                          pDataForLabel = NULL,
                          pDataForShape = "type",
                          seed = 0)
    
    vdiffr::expect_doppelganger(
        "ggplotSamplesMDS with axes 1 and 2 and extVars",
        fig = p)
    
    p <- ggplotSamplesMDS(mdsObj = mdsObj,
                          pData = flowCore::pData(fsAll),
                          projectionAxes = c(3,4),
                          biplot = TRUE,
                          extVariables = extVars,
                          pDataForColour = "group",
                          pDataForLabel = "name",
                          seed = 0)
    
    vdiffr::expect_doppelganger(
        "ggplotSamplesMDS with axes 3 and 4 and extVars",
        fig = p)
    
    p <- ggplotSamplesMDS(mdsObj = mdsObj,
                          pData = flowCore::pData(fsAll),
                          projectionAxes = c(3,4),
                          biplot = TRUE,
                          extVariables = extVars,
                          pDataForColour = "group",
                          pDataForLabel = NULL,
                          seed = 0)
    
    vdiffr::expect_doppelganger(
        "ggplotSamplesMDS no point labels",
        fig = p)
    
})


test_that("ggplotSamplesMDSShepard works", {
    transList <- estimateScaleTransforms(
        ff = OMIP021Samples[[1]],
        fluoMethod = "estimateLogicle",
        scatterMethod = "linearQuantile",
        scatterRefMarker = "BV785 - CD3")
    
    ffList <- flowCore::flowSet_to_list(OMIP021Samples)
    
    for(i in 3:5){
        ffList[[i]] <- 
            aggregateAndSample(
                OMIP021Samples,
                seed = 10*i,
                nTotalEvents = 5000)[,1:22]
    }
    
    fsNames <- c("Donor1", "Donor2", paste0("Agg",1:3))
    names(ffList) <- fsNames
    
    fsAll <- as(ffList,"flowSet")
    
    flowCore::pData(fsAll)$group <- factor(c("D1", "D2", rep("Agg", 3)))
    flowCore::pData(fsAll)$type <- factor(c("real", "real", rep("synthetic", 3)))
    
    pwDist <- getPairWiseEMDDist(fsAll, 
                                 channels = c("FSC-A", "SSC-A"),
                                 transList = transList,
                                 verbose = FALSE)
    
    mdsObj <- computeMetricMDS(pwDist, nDim = 4, seed = 0)
    
    p <- ggplotSamplesMDSShepard(mdsObj,
                                 nDim = 2,
                                 pointSize = 1,
                                 title = "Shepard with 2 dimensions")
    
    vdiffr::expect_doppelganger(
        "ggplotSamplesMDSShepard with 2 dimensions",
        fig = p)
    
    p <- ggplotSamplesMDSShepard(mdsObj,
                                 nDim = 3,
                                 title = "Shepard with 3 dimensions") 
    
    vdiffr::expect_doppelganger(
        "ggplotSamplesMDSShepard with 3 dimensions",
        fig = p)
    
    p <- ggplotSamplesMDSShepard(mdsObj,
                                 title = "Shepard with default nb of dimensions") 
    
    vdiffr::expect_doppelganger(
        "ggplotSamplesMDSShepard with default dim nb",
        fig = p)
    
})

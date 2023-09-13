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

test_that("getEMDDist works", {
    transList <- estimateScaleTransforms(
        ff = OMIP021Samples[[1]],
        fluoMethod = "estimateLogicle",
        scatterMethod = "linearQuantile",
        scatterRefMarker = "BV785 - CD3")
    
    # distance with itself, all channels
    distDum <- getEMDDist(ff1 = OMIP021Samples[[1]],
                          ff2 = OMIP021Samples[[1]],
                          transList = transList,
                          binSize = 0.05,
                          minRange = -10,
                          maxRange = 10,
                          returnAll = FALSE)
    expect_equal(distDum, 0.)
    
    # returning only distance, 2 channels, with transList
    dist1 <- getEMDDist(ff1 = OMIP021Samples[[1]], 
                        ff2 = OMIP021Samples[[2]], 
                        channels = c("FSC-A", "SSC-A"),
                        transList = transList,
                        binSize = 0.05,
                        minRange = -10,
                        maxRange = 10,
                        returnAll = FALSE)
    
    expect_equal(dist1, 0.1551)
    
    
    # returning only distance, 2 channels, no transList
    ffTrans1 <- flowCore::transform(OMIP021Samples[[1]], transList)
    ffTrans2 <- flowCore::transform(OMIP021Samples[[2]], transList)
    dist2 <- getEMDDist(ff1 = ffTrans1, 
                        ff2 = ffTrans2, 
                        channels = c("FSC-A", "SSC-A"),
                        binSize = 0.05,
                        minRange = -10,
                        maxRange = 10,
                        returnAll = FALSE)
    
    expect_equal(dist2, dist1)
    
    # using only one channel, passed by marker name
    dist3 <- getEMDDist(ff1 = OMIP021Samples[[1]], 
                        ff2 = OMIP021Samples[[2]], 
                        channels = c("BV785 - CD3"),
                        transList = transList,
                        binSize = 0.05,
                        minRange = -10,
                        maxRange = 10,
                        returnAll = FALSE)
    
    expect_equal(dist3, 0.1393)
    
    # using only one channel, passed by index
    dist4 <- getEMDDist(ff1 = OMIP021Samples[[1]], 
                        ff2 = OMIP021Samples[[2]], 
                        channels = 10,
                        transList = transList,
                        binSize = 0.05,
                        minRange = -10,
                        maxRange = 10,
                        returnAll = FALSE)
    
    expect_equal(dist4, dist3)
    
    # check that a warning is issued, when [minRange, maxRange] does not span
    # all events
    w <- capture_warnings({
        distWarn <- getEMDDist(ff1 = OMIP021Samples[[1]],
                               ff2 = OMIP021Samples[[1]],
                               channels = c("FSC-A", "SSC-A"),
                               transList = transList,
                               binSize = 0.05,
                               minRange = -3,
                               maxRange = 3,
                               returnAll = FALSE)
    })
    expect_equal(length(w), 4)
    for(i in 1:4){
        expect_match(w[i], regexp = "does not span all events for channel")  
    }
    expect_equal(distWarn, 0.)
    
    #returning all
    allDist <- getEMDDist(ff1 = OMIP021Samples[[1]], 
                          ff2 = OMIP021Samples[[2]], 
                          channels = c("FSC-A", "SSC-A"),
                          transList = transList,
                          binSize = 0.05,
                          minRange = -10,
                          maxRange = 10,
                          returnAll = TRUE)
    expectedDists <- c(0.11292, 0.04218)
    names(expectedDists) <- c("FSC-A", "SSC-A")
    expect_equal(allDist$distances, expectedDists)
    
    
})

test_that("getPairWiseEMDDist works", {
    transList <- CytoPipeline::estimateScaleTransforms(
        ff = OMIP021Samples[[1]],
        fluoMethod = "estimateLogicle",
        scatterMethod = "linearQuantile",
        scatterRefMarker = "BV785 - CD3")
    
    pwDist <- getPairWiseEMDDist(fs = OMIP021Samples,
                                 channels = c("FSC-A", "SSC-A"),
                                 transList = transList,
                                 binSize = 0.05,
                                 minRange = -10,
                                 maxRange = 10
    )
    expect_equal(dim(pwDist), c(2,2))
    expect_equal(pwDist[1,1], 0.)
    expect_equal(pwDist[1,2], 0.1551)
    expect_equal(pwDist[2,1], 0.1551)
    expect_equal(pwDist[2,2], 0.)
    
})

test_that("getChannelsSummaryStat works", {
    
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
    
    channelsOrMarkers <- c("FSC-A", "SSC-A", "BV785 - CD3")
    
    ret <- getChannelsSummaryStat(
        fsAll,
        channels = channelsOrMarkers,
        transList = transList,
        statFUN = stats::median,
        verbose = FALSE)
    
    expect_equal(unname(colnames(ret)), channelsOrMarkers)
    
    expect_equal(unname(ret[2,1]), 1.941572819633120)
    expect_equal(unname(ret[3,2]), 1.131832378296435)
    expect_equal(unname(ret[4,3]), 1.638696159955064)
    
    ret <- getChannelsSummaryStat(
        fsAll,
        channels = channelsOrMarkers,
        transList = transList,
        statFUN = mean,
        verbose = FALSE)
    
    expect_equal(unname(colnames(ret)), channelsOrMarkers)
    
    expect_equal(unname(ret[2,1]), 1.961771149533659)
    expect_equal(unname(ret[3,2]), 1.393034241892733)
    expect_equal(unname(ret[4,3]), 1.907561024702794)
    
    ret <- getChannelsSummaryStat(
        fsAll,
        channels = channelsOrMarkers,
        transList = transList,
        statFUN = stats::sd,
        verbose = FALSE)
    
    expect_equal(unname(colnames(ret)), channelsOrMarkers)
    
    expect_equal(unname(ret[2,1]), 0.5942291770870205)
    expect_equal(unname(ret[3,2]), 0.7352746696905406)
    expect_equal(unname(ret[4,3]), 0.8420740225208073)
    
})

test_that("computeMetricMDS works", {
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
    
    pwDist <- getPairWiseEMDDist(fsAll, 
                                 channels = c("FSC-A", "SSC-A"),
                                 transList = transList,
                                 verbose = FALSE)
    
    mdsObj <- computeMetricMDS(pwDist, seed = 0)
    
    expect_equal(mdsObj$stress, 0.0203635387)
    tgtspp <- c(
        11.111867, 13.086625, 4.448069, 37.334795, 34.018645)
    names(tgtspp) <- 1:5
    expect_equal(mdsObj$spp, tgtspp)
    
})

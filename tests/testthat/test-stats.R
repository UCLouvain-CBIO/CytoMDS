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

outputDir <- base::tempdir()

transList <- estimateScaleTransforms(
    ff = OMIP021Samples[[1]],
    fluoMethod = "estimateLogicle",
    scatterMethod = "linearQuantile",
    scatterRefMarker = "BV785 - CD3")

OMIP021Trans <- CytoPipeline::applyScaleTransforms(
    OMIP021Samples, 
    transList)

test_that("getEMDDist works", {
    # distance with itself, all channels
    distDum <- getEMDDist(ff1 = OMIP021Trans[[1]],
                          ff2 = OMIP021Trans[[1]],
                          binSize = 0.05,
                          minRange = -10,
                          maxRange = 10,
                          returnAll = FALSE)
    expect_equal(distDum, 0.)
    
    # returning only distance, 2 channels
    dist1 <- getEMDDist(ff1 = OMIP021Trans[[1]], 
                        ff2 = OMIP021Trans[[2]], 
                        channels = c("FSC-A", "SSC-A"),
                        binSize = 0.05,
                        minRange = -10,
                        maxRange = 10,
                        returnAll = FALSE)
    
    expect_equal(dist1, 0.1551)
    
    # using only one channel, passed by marker name
    dist3 <- getEMDDist(ff1 = OMIP021Trans[[1]], 
                        ff2 = OMIP021Trans[[2]], 
                        channels = c("BV785 - CD3"),
                        binSize = 0.05,
                        minRange = -10,
                        maxRange = 10,
                        returnAll = FALSE)
    
    expect_equal(dist3, 0.1393)
    
    # using only one channel, passed by index
    dist4 <- getEMDDist(ff1 = OMIP021Trans[[1]], 
                        ff2 = OMIP021Trans[[2]], 
                        channels = 10,
                        binSize = 0.05,
                        minRange = -10,
                        maxRange = 10,
                        returnAll = FALSE)
    
    expect_equal(dist4, dist3)
    
    # check that a warning is issued, when [minRange, maxRange] does not span
    # all events
    w <- capture_warnings({
        distWarn <- getEMDDist(ff1 = OMIP021Trans[[1]],
                               ff2 = OMIP021Trans[[1]],
                               channels = c("FSC-A", "SSC-A"),
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
    allDist <- getEMDDist(ff1 = OMIP021Trans[[1]], 
                          ff2 = OMIP021Trans[[2]], 
                          channels = c("FSC-A", "SSC-A"),
                          binSize = 0.05,
                          minRange = -10,
                          maxRange = 10,
                          returnAll = TRUE)
    expectedDists <- c(0.11292, 0.04218)
    names(expectedDists) <- c("FSC-A", "SSC-A")
    expect_equal(allDist$distances, expectedDists)
    
    
})

test_that("getPairwiseEMDDist works", {
    pwDist <- getPairwiseEMDDist(fs = OMIP021Trans,
                                 channels = c("FSC-A", "SSC-A"),
                                 binSize = 0.05,
                                 minRange = -10,
                                 maxRange = 10
    )
    expect_equal(dim(pwDist), c(2,2))
    expect_equal(pwDist[1,1], 0.)
    expect_equal(pwDist[1,2], 0.1551)
    expect_equal(pwDist[2,1], 0.1551)
    expect_equal(pwDist[2,2], 0.)
    
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
    
    pwDist2 <- getPairwiseEMDDist(
        fs = fsAll[1:3],
        fs2 = fsAll[4:5],
        channels = c("FSC-A", "SSC-A"),
        binSize = 0.05,
        minRange = -10,
        maxRange = 10)
    expect_equal(dim(pwDist2), c(3,2))
    expect_equal(pwDist2[1,1], 0.07241)
    expect_equal(pwDist2[1,2], 0.08347)
    expect_equal(pwDist2[2,1], 0.08501)
    expect_equal(pwDist2[2,2], 0.07451)
    expect_equal(pwDist2[3,1], 0.01293)
    expect_equal(pwDist2[3,2], 0.01813)
})

test_that("getPairwiseEMDDist works with BiocParallel", {
    
    logDir <- file.path(outputDir, "BiocParallel", "log")
    
    suppressWarnings(dir.create(logDir, recursive = TRUE))
    
    bp <- BiocParallel::SnowParam(workers = 2, log = TRUE, 
                                  logdir = logDir,
                                  progressbar = TRUE)
    
    pwDist <- suppressWarnings(getPairwiseEMDDist(
        fs = OMIP021Trans,
        channels = c("FSC-A", "SSC-A"),
        binSize = 0.05,
        minRange = -10,
        maxRange = 10,
        useBiocParallel = TRUE,
        BPPARAM = bp))
    
    expect_equal(dim(pwDist), c(2,2))
    expect_equal(pwDist[1,1], 0.)
    expect_equal(pwDist[1,2], 0.1551)
    expect_equal(pwDist[2,1], 0.1551)
    expect_equal(pwDist[2,2], 0.)
    
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
    
    expect_error (
        getPairwiseEMDDist(
            fs = fsAll[1:3],
            fs2 = fsAll[4:5],
            channels = c("FSC-A", "SSC-A"),
            binSize = 0.05,
            minRange = -10,
            maxRange = 10,
            useBiocParallel = TRUE,
            BPPARAM = bp), 
        regexp = "Using BiocParallel currently only works")
    
    pwDist2 <- suppressWarnings(getPairwiseEMDDist(
        fs = fsAll[1:2],
        fs2 = fsAll[4:5],
        channels = c("FSC-A", "SSC-A"),
        binSize = 0.05,
        minRange = -10,
        maxRange = 10,
        useBiocParallel = TRUE,
        BPPARAM = bp))
        
    expect_equal(dim(pwDist2), c(2,2))
    expect_equal(pwDist2[1,1], 0.07241)
    expect_equal(pwDist2[1,2], 0.08347)
    expect_equal(pwDist2[2,1], 0.08501)
    expect_equal(pwDist2[2,2], 0.07451)
})

test_that("getChannelSummaryStats works", {
   
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
    
    channelsOrMarkers <- c("FSC-A", "SSC-A", "BV785 - CD3")
    
    ret <- getChannelSummaryStats(
        fsAll,
        channels = channelsOrMarkers,
        statFUNs = stats::median,
        verbose = FALSE)
    
    expect_equal(unname(colnames(ret)), channelsOrMarkers)
    
    expect_equal(unname(ret[2,1]), 1.941572819633120)
    expect_equal(unname(ret[3,2]), 1.131832378296435)
    expect_equal(unname(ret[4,3]), 1.638696159955064)
    
    ret <- getChannelSummaryStats(
        fsAll,
        channels = channelsOrMarkers,
        statFUNs = list("mean" = mean, "std.dev" = stats::sd),
        verbose = FALSE)
    
    expect_equal(names(ret), c("mean", "std.dev"))
    expect_equal(unname(colnames(ret[[1]])), channelsOrMarkers)
    expect_equal(unname(ret[[1]][2,1]), 1.961771149533659)
    expect_equal(unname(ret[[1]][3,2]), 1.393034241892733)
    expect_equal(unname(ret[[1]][4,3]), 1.907561024702794)
    
    
    expect_equal(unname(colnames(ret[[2]])), channelsOrMarkers)
    
    expect_equal(unname(ret[[2]][2,1]), 0.5942291770870205)
    expect_equal(unname(ret[[2]][3,2]), 0.7352746696905406)
    expect_equal(unname(ret[[2]][4,3]), 0.8420740225208073)
    
})

test_that("computeMetricMDS works", {
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
    
    pwDist <- getPairwiseEMDDist(fsAll, 
                                 channels = c("FSC-A", "SSC-A"),
                                 verbose = FALSE)
    
    mdsObj <- computeMetricMDS(pwDist, nDim = 2, seed = 0)
    
    expect_equal(mdsObj$stress, 0.0203635387)
    tgtspp <- c(
        11.111867, 13.086625, 4.448069, 37.334795, 34.018645)
    names(tgtspp) <- 1:5
    expect_equal(mdsObj$spp, tgtspp)
    
    # with no user provided nDim, but (implicit) target pseudo rsquare = 0.95
    mdsObj2 <- computeMetricMDS(pwDist, seed = 0)
    
    expect_equal(mdsObj$nDim, 2)
    expect_equal(mdsObj$RSq[2], 0.99843722)
    
    # with no user provided nDim, but explicit target pseudo rsquare = 0.98
    mdsObj3 <- computeMetricMDS(pwDist, seed = 0, targetPseudoRSq = 0.98)
    
    expect_equal(mdsObj3$nDim, 2)
    expect_equal(mdsObj3$RSq[2], 0.99843722)
    
    # with no user provided nDim, but explicit target pseudo rsquare = 0.999
    mdsObj4 <- computeMetricMDS(pwDist, seed = 0, targetPseudoRSq = 0.999)
    
    expect_equal(mdsObj4$nDim, 3)
    expect_equal(mdsObj4$RSq[3], 0.99988906)
    
})

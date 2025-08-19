# CytoMDS - Copyright (C) <2023-2025>
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

# function used in subsequent tests where the flow set is not stored at once
# in memory but the flow frames are loaded dynamically upon request
simulMemoryLoad <- function(
        exprMatrixIndex, 
        filePaths, 
        nSamples,
        channels = NULL,
        verbose) {
    
    if (!is.numeric(nSamples) || nSamples < 1) {
        stop("nSamples should be a numeric greater than or equal to 1")
    }
    
    if (!is.numeric(exprMatrixIndex) || exprMatrixIndex < 1) {
        stop("exprMatrixIndex should be a numeric greater than or equal to 1")
    }
    
    if (exprMatrixIndex > nSamples) {
        stop("exprMatrixIndex should be <= nSamples")
    }
    
    # init Flow Set
    if (verbose) {
        message("loading ff #", exprMatrixIndex, "/", nSamples, "...")
        message("reading initial flow set...")
    }
    
    initialFS <- flowCore::read.flowSet(
        filePaths,
        transformation = "linearize",
        alter.names = FALSE,
        name = "OMIP-21",
        truncate_max_range = TRUE,
        min.limit = NULL)
    
    if (verbose) {
        message("transforming samples...")
    }
    
    transList <- estimateScaleTransforms(
        ff = initialFS[[1]],
        fluoMethod = "estimateLogicle",
        scatterMethod = "linearQuantile",
        scatterRefMarker = "BV785 - CD3")
    
    transFS <- CytoPipeline::applyScaleTransforms(
        initialFS, 
        transList)
    
    if (verbose) {
        message("creating specific sample...")
    }
    
    if (exprMatrixIndex <= ceiling(nSamples/2)){
        ff <- CytoPipeline::subsample(
            transFS[[1]],
            nEvents = 1000,
            seed = exprMatrixIndex)[,1:22]
    } else {
        ff <- CytoPipeline::subsample(
            transFS[[2]],
            nEvents = 1000,
            seed = exprMatrixIndex)[,1:22]
    }
    
    # fscMean <- mean(flowCore::exprs(ff)[,"FSC-A"])
    # sscMean <- mean(flowCore::exprs(ff)[,"SSC-A"])
    # message("OBTAINED FF: FSC mean: ", round(fscMean,6), "; SSC mean:", 
    #          round(sscMean,6), "\n")
    
    exprMat <- flowCore::exprs(ff)
    channelsId <- CytoMDS:::.toChannelNames(channels, ff)
    
    exprMat <- exprMat[, channelsId, drop=FALSE]
    if (!is.null(channels)) {
        colnames(exprMat) <- channels
    }
    exprMat
}

test_that("EMDDist works with ffs", {
    # distance with itself, all channels
    distDum <- EMDDist(x1 = OMIP021Trans[[1]],
                       x2 = OMIP021Trans[[1]],
                       binSize = 0.05,
                       minRange = -10,
                       maxRange = 10,
                       returnAll = FALSE)
    expect_equal(distDum, 0.)
    
    # returning only distance, 2 channels
    dist1 <- EMDDist(x1 = OMIP021Trans[[1]], 
                     x2 = OMIP021Trans[[2]], 
                     channels = c("FSC-A", "SSC-A"),
                     binSize = 0.05,
                     minRange = -10,
                     maxRange = 10,
                     returnAll = FALSE)
    
    expect_equal(dist1, 0.1551)
    
    # using only one channel, passed by marker name
    dist3 <- EMDDist(x1 = OMIP021Trans[[1]], 
                     x2 = OMIP021Trans[[2]], 
                     channels = c("BV785 - CD3"),
                     binSize = 0.05,
                     minRange = -10,
                     maxRange = 10,
                     returnAll = FALSE)
    
    expect_equal(dist3, 0.1393)
    
    # using only one channel, passed by index
    dist4 <- EMDDist(x1 = OMIP021Trans[[1]], 
                     x2 = OMIP021Trans[[2]], 
                     channels = 10,
                     binSize = 0.05,
                     minRange = -10,
                     maxRange = 10,
                     returnAll = FALSE)
    
    expect_equal(dist4, dist3)
    
    # check that a warning is issued, when [minRange, maxRange] does not span
    # all events
    w <- capture_warnings({
        distWarn <- EMDDist(x1 = OMIP021Trans[[1]],
                            x2 = OMIP021Trans[[1]],
                            channels = c("FSC-A", "SSC-A"),
                            binSize = 0.05,
                            minRange = -3,
                            maxRange = 3,
                            returnAll = FALSE)
    })
    expect_equal(length(w), 4)
    expect_match(w[1], regexp = "does not span all events for channel")  
    expect_equal(distWarn, 0.)
    
    #returning all
    allDist <- EMDDist(x1 = OMIP021Trans[[1]], 
                       x2 = OMIP021Trans[[2]], 
                       channels = c("FSC-A", "SSC-A"),
                       binSize = 0.05,
                       minRange = -10,
                       maxRange = 10,
                       returnAll = TRUE)
    expectedDists <- c(0.11292, 0.04218)
    names(expectedDists) <- c("FSC-A", "SSC-A")
    expect_equal(allDist$distances, expectedDists)
})

test_that("EMDDist works with expr matrices", {
    expr1 <- flowCore::exprs(OMIP021Trans[[1]])[, 
        CytoPipeline::areSignalCols(OMIP021Trans[[1]])]
        
    expr2 <- flowCore::exprs(OMIP021Trans[[2]])[, 
        CytoPipeline::areSignalCols(OMIP021Trans[[2]])]
    
    # distance with itself, all channels
    distDum <- EMDDist(x1 = expr1,
                       x2 = expr1,
                       binSize = 0.05,
                       minRange = -10,
                       maxRange = 10,
                       returnAll = FALSE)
    expect_equal(distDum, 0.)
    
    # returning only distance, 2 channels
    dist1 <- EMDDist(x1 = expr1, 
                     x2 = expr2, 
                     channels = c("FSC-A", "SSC-A"),
                     binSize = 0.05,
                     minRange = -10,
                     maxRange = 10,
                     returnAll = FALSE)
    
    expect_equal(dist1, 0.1551)
    
    # using only one channel, passed by index
    dist4 <- EMDDist(x1 = expr1, 
                     x2 = expr2, 
                     channels = 10,
                     binSize = 0.05,
                     minRange = -10,
                     maxRange = 10,
                     returnAll = FALSE)
    
    expect_equal(dist4, 0.1393)
    
    # check that a warning is issued, when [minRange, maxRange] does not span
    # all events
    w <- capture_warnings({
        distWarn <- EMDDist(x1 = expr1,
                            x2 = expr1,
                            channels = c("FSC-A", "SSC-A"),
                            binSize = 0.05,
                            minRange = -3,
                            maxRange = 3,
                            returnAll = FALSE)
    })
    expect_equal(length(w), 4)
    expect_match(w[1], regexp = "does not span all events for channel")  
    expect_equal(distWarn, 0.)
    
    #returning all
    allDist <- EMDDist(x1 = expr1, 
                       x2 = expr2, 
                       channels = c("FSC-A", "SSC-A"),
                       binSize = 0.05,
                       minRange = -10,
                       maxRange = 10,
                       returnAll = TRUE)
    expectedDists <- c(0.11292, 0.04218)
    names(expectedDists) <- c("FSC-A", "SSC-A")
    expect_equal(allDist$distances, expectedDists)
})


test_that("generateBlocks2D works", {
    
    # trivial cases with 1 sample
    blocks2D <- .generateBlocks2D(
        rowRange = c(1, 1),
        colRange = c(1, 1),
        nRowBlock = 1)
    
    expect_equal(length(blocks2D), 0)
    # expect_equal(blocks2D[[1]]$rowMin, 1)
    # expect_equal(blocks2D[[1]]$rowMax, 1)
    # expect_equal(blocks2D[[1]]$colMin, 1)
    # expect_equal(blocks2D[[1]]$colMax, 1)
    
    # cases that imposes the nb of blocks in 1D
    blocks2D <- .generateBlocks2D(
        rowRange = c(1,10),
        colRange = c(1,10),
        nRowBlock = 2)
    
    expect_equal(length(blocks2D), 3)
    expect_equal(blocks2D[[2]]$rowMin, 1)
    expect_equal(blocks2D[[2]]$rowMax, 5)
    expect_equal(blocks2D[[2]]$colMin, 6)
    expect_equal(blocks2D[[2]]$colMax, 10)
    
    blocks2D <- .generateBlocks2D(
        rowRange = c(1,10),
        colRange = c(1,10),
        nRowBlock = 3)
    
    expect_equal(length(blocks2D), 6)
    expect_equal(blocks2D[[2]]$rowMin, 1)
    expect_equal(blocks2D[[2]]$rowMax, 4)
    expect_equal(blocks2D[[2]]$colMin, 5)
    expect_equal(blocks2D[[2]]$colMax, 7)
    
    blocks2D <- .generateBlocks2D(
        rowRange = c(1,10),
        colRange = c(9,10),
        nRowBlock = 3)
    
    expect_equal(length(blocks2D),3)
    expect_equal(blocks2D[[2]]$rowMin, 5)
    expect_equal(blocks2D[[2]]$rowMax, 7)
    expect_equal(blocks2D[[2]]$colMin, 9)
    expect_equal(blocks2D[[2]]$colMax, 10)
    
    blocks2D <- .generateBlocks2D(
        rowRange = c(1,10),
        colRange = c(1,22),
        nRowBlock = 1)
    
    expect_equal(length(blocks2D),3)
    expect_equal(blocks2D[[2]]$rowMin, 1)
    expect_equal(blocks2D[[2]]$rowMax, 10)
    expect_equal(blocks2D[[2]]$colMin, 9)
    expect_equal(blocks2D[[2]]$colMax, 15)
    
    blocks2D <- .generateBlocks2D(
        rowRange = c(1,10),
        colRange = c(10,10),
        nRowBlock = 10)
    
    expect_equal(length(blocks2D),9)
    expect_equal(blocks2D[[2]]$rowMin, 2)
    expect_equal(blocks2D[[2]]$rowMax, 2)
    expect_equal(blocks2D[[2]]$colMin, 10)
    expect_equal(blocks2D[[2]]$colMax, 10)
    
    blocks2D <- .generateBlocks2D(
        rowRange = c(1,10),
        colRange = c(11,20),
        nRowBlock = 2)
    
    expect_equal(length(blocks2D),4)
    expect_equal(blocks2D[[2]]$rowMin, 1)
    expect_equal(blocks2D[[2]]$rowMax, 5)
    expect_equal(blocks2D[[2]]$colMin, 16)
    expect_equal(blocks2D[[2]]$colMax, 20)
    
    blocks2D <- .generateBlocks2D(
        rowRange = c(1,10),
        colRange = c(11,21),
        nRowBlock = 2)
    
    expect_equal(length(blocks2D),6)
    expect_equal(blocks2D[[2]]$rowMin, 1)
    expect_equal(blocks2D[[2]]$rowMax, 5)
    expect_equal(blocks2D[[2]]$colMin, 15)
    expect_equal(blocks2D[[2]]$colMax, 17)
    
    blocks2D <- .generateBlocks2D(
        rowRange = c(1,1),
        colRange = c(1,10),
        nRowBlock = 1)
    
    expect_equal(length(blocks2D),9)
    expect_equal(blocks2D[[2]]$rowMin, 1)
    expect_equal(blocks2D[[2]]$rowMax, 1)
    expect_equal(blocks2D[[2]]$colMin, 3)
    expect_equal(blocks2D[[2]]$colMax, 3)
    
    blocks2D <- .generateBlocks2D(
        rowRange = c(1,2), 
        colRange = c(1,2), 
        nRowBlock = 2)
    
    expect_equal(length(blocks2D),1)
    expect_equal(blocks2D[[1]]$rowMin, 1)
    expect_equal(blocks2D[[1]]$rowMax, 1)
    expect_equal(blocks2D[[1]]$colMin, 2)
    expect_equal(blocks2D[[1]]$colMax, 2)
    
})

test_that("optimizeRowBlockSize works", {
    
    nRowBlock <- .optimizeRowBlockNb(
        rowRange = c(1,1), 
        colRange = c(1,10), 
        nCores = 1, 
        memSize = 1)
    
    expect_equal(nRowBlock, 1)
    
    nRowBlock <- .optimizeRowBlockNb(
        rowRange = c(1,10), 
        colRange = c(10,10), 
        nCores = 1, 
        memSize = 1)
    
    expect_equal(nRowBlock, 10)
    
    nRowBlock <- .optimizeRowBlockNb(
        rowRange = c(1,10), 
        colRange = c(10,10), 
        nCores = 1)
    
    expect_equal(nRowBlock, 1)
    
    nRowBlock <- .optimizeRowBlockNb(
        rowRange = c(1,10), 
        colRange = c(10,10), 
        nCores = 1,
        memSize = 10)
    
    expect_equal(nRowBlock, 1)
    
    nRowBlock <- .optimizeRowBlockNb(
        rowRange = c(1,10), 
        colRange = c(10,10), 
        nCores = 1,
        memSize = 6)
    
    expect_equal(nRowBlock, 2)
    
    nRowBlock <- .optimizeRowBlockNb(
        rowRange = c(1,10), 
        colRange = c(10,10), 
        nCores = 4,
        memSize = 6)
    
    expect_equal(nRowBlock, 4)
    
    nRowBlock <- .optimizeRowBlockNb(
        rowRange = c(1,10), 
        colRange = c(10,10), 
        nCores = 4,
        memSize = 2)
    
    expect_equal(nRowBlock, 10)
    
    nRowBlock <- .optimizeRowBlockNb(
        rowRange = c(1,12),
        colRange = c(1,12),
        nCores = 8,
        memSize = 6 
    )
    
    expect_equal(nRowBlock, 4)
    
    nRowBlock <- .optimizeRowBlockNb(
        rowRange = c(1,12),
        colRange = c(1,12),
        nCores = 8,
        memSize = 12
    )
    
    expect_equal(nRowBlock, 4)
    
    nRowBlock <- .optimizeRowBlockNb(
        rowRange = c(1,12),
        colRange = c(1,12),
        memSize = 12
    )
    
    expect_equal(nRowBlock, 1)
    
    })

test_that("pairwiseEMDDist with fs works", {
    pwDist <- pairwiseEMDDist(
        x = OMIP021Trans[1],
        channels = c("FSC-A", "SSC-A"),
        binSize = 0.05,
        minRange = -10,
        maxRange = 10)
    expect_equal(dim(pwDist), c(1,1))
    expect_equal(as.matrix(pwDist)[1,1], 0.)
    
    pwDist <- pairwiseEMDDist(
        x = OMIP021Trans,
        channels = c("FSC-A", "SSC-A"),
        binSize = 0.05,
        minRange = -10,
        maxRange = 10)
    
    pwDist <- as.matrix(pwDist)
    
    expect_equal(dim(pwDist), c(2,2))
    expect_equal(pwDist[1,1], 0.)
    expect_equal(pwDist[1,2], 0.1551)
    expect_equal(pwDist[2,1], 0.1551)
    expect_equal(pwDist[2,2], 0.)
    
    # now calculate the distance with all channels
    # and select later
    pwDist <- pairwiseEMDDist(
        x = OMIP021Trans,
        channels = NULL,
        binSize = 0.05,
        minRange = -10,
        maxRange = 10)
    
    pwDist <- as.matrix(pwDist, whichFeatures = c("FSC-A", "SSC-A"))
    expect_equal(dim(pwDist), c(2,2))
    expect_equal(pwDist[1,1], 0.)
    expect_equal(pwDist[1,2], 0.1551)
    expect_equal(pwDist[2,1], 0.1551)
    expect_equal(pwDist[2,2], 0.)
    
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
    
    pwDist2 <- pairwiseEMDDist(
        x = fsAll,
        rowRange = c(1,3),
        colRange = c(4,5),
        channels = c("FSC-A", "SSC-A"),
        binSize = 0.05,
        minRange = -10,
        maxRange = 10,
        verbose = FALSE)
    
    pwDist2 <- as.matrix(pwDist2)
    
    expect_equal(dim(pwDist2), c(3,2))
    expect_equal(pwDist2[1,1], 0.07241)
    expect_equal(pwDist2[1,2], 0.08347)
    expect_equal(pwDist2[2,1], 0.08501)
    expect_equal(pwDist2[2,2], 0.07451)
    expect_equal(pwDist2[3,1], 0.01293)
    expect_equal(pwDist2[3,2], 0.01813)
    
    # test with asymetric block overlapping with diagonal
    pwDist3 <- pairwiseEMDDist(
        x = fsAll,
        rowRange = c(1,5),
        colRange = c(2,5),
        channels = c("FSC-A", "SSC-A"),
        binSize = 0.05,
        minRange = -10,
        maxRange = 10,
        verbose = FALSE)
    
    pwDist3 <- as.matrix(pwDist3)
    
    expect_equal(dim(pwDist3), c(5,4))
    expect_equal(pwDist3[1,1], 0.1551)
    expect_equal(pwDist3[1,2], 0.07132)
    expect_equal(pwDist3[2,1], 0.)
    expect_equal(pwDist3[2,2], 0.08556)
    expect_equal(pwDist3[3,1], 0.)
    expect_equal(pwDist3[3,2], 0.)
    
    # with verbose
    msg <- capture_messages(pwDist4 <- pairwiseEMDDist(
        x = fsAll,
        rowRange = c(1,5),
        colRange = c(2,5),
        channels = c("FSC-A", "SSC-A"),
        binSize = 0.05,
        minRange = -10,
        maxRange = 10,
        verbose = TRUE))
    
    expect_equal(length(msg), 22)
    expect_equal(msg[5], "Calculating histogram for file 2...\n")
    expect_equal(msg[10], "Loading file 5...\n")
    expect_equal(msg[21], "i = 3; j = 5; sum(dist) = 0.01813\n")
    
    pwDist4 <- as.matrix(pwDist4)
    
    expect_equal(dim(pwDist4), c(5,4))
    expect_equal(pwDist4[1,1], 0.1551)
    expect_equal(pwDist4[1,2], 0.07132)
    expect_equal(pwDist4[2,1], 0.)
    expect_equal(pwDist4[2,2], 0.08556)
    expect_equal(pwDist4[3,1], 0.)
    expect_equal(pwDist4[3,2], 0.)
})

test_that("pairwiseEMDDist with expr matrix works", {
    
    expr1 <- flowCore::exprs(OMIP021Trans[[1]])[, 
        CytoPipeline::areSignalCols(OMIP021Trans[[1]])]
    expr2 <- flowCore::exprs(OMIP021Trans[[2]])[, 
        CytoPipeline::areSignalCols(OMIP021Trans[[2]])]
    
    pwDist <- pairwiseEMDDist(
        x = list(expr1),
        channels = c("FSC-A", "SSC-A"),
        binSize = 0.05,
        minRange = -10,
        maxRange = 10)
    
    pwDist <- as.matrix(pwDist)
    
    expect_equal(dim(pwDist), c(1,1))
    expect_equal(pwDist[1,1], 0.)
    
    exprList <- list(sample1 = expr1, sample2 = expr2)
    
    pwDist <- pairwiseEMDDist(
        x = exprList,
        channels = c("FSC-A", "SSC-A"),
        binSize = 0.05,
        minRange = -10,
        maxRange = 10)
    
    pwDist <- as.matrix(pwDist)
    
    expect_equal(dim(pwDist), c(2,2))
    expect_equal(pwDist[1,1], 0.)
    expect_equal(pwDist[1,2], 0.1551)
    expect_equal(pwDist[2,1], 0.1551)
    expect_equal(pwDist[2,2], 0.)
    
    # now calculate the distance with all channels
    # and select later
    pwDist <- pairwiseEMDDist(
        x = list(expr1, expr2),
        channels = NULL,
        binSize = 0.05,
        minRange = -10,
        maxRange = 10)
    pwDist <- as.matrix(pwDist, whichFeatures = c("FSC-A", "SSC-A"))
    expect_equal(dim(pwDist), c(2,2))
    expect_equal(pwDist[1,1], 0.)
    expect_equal(pwDist[1,2], 0.1551)
    expect_equal(pwDist[2,1], 0.1551)
    expect_equal(pwDist[2,2], 0.)
    
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
    
    exprList <- lapply(ffList, FUN = function(ff) {
        expr <- flowCore::exprs(ff)[, CytoPipeline::areSignalCols(ff)]
    })
        
    pwDist2 <- pairwiseEMDDist(
        x = exprList,
        rowRange = c(1,3),
        colRange = c(4,5),
        channels = c("FSC-A", "SSC-A"),
        binSize = 0.05,
        minRange = -10,
        maxRange = 10,
        verbose = FALSE)
    
    pwDist2 <- as.matrix(pwDist2)
    
    expect_equal(dim(pwDist2), c(3,2))
    expect_equal(pwDist2[1,1], 0.07241)
    expect_equal(pwDist2[1,2], 0.08347)
    expect_equal(pwDist2[2,1], 0.08501)
    expect_equal(pwDist2[2,2], 0.07451)
    expect_equal(pwDist2[3,1], 0.01293)
    expect_equal(pwDist2[3,2], 0.01813)
    
    # test with asymetric block overlapping with diagonal
    pwDist3 <- pairwiseEMDDist(
        x = exprList,
        rowRange = c(1,5),
        colRange = c(2,5),
        channels = c("FSC-A", "SSC-A"),
        binSize = 0.05,
        minRange = -10,
        maxRange = 10,
        verbose = FALSE)
    
    pwDist3 <- as.matrix(pwDist3)
    
    expect_equal(dim(pwDist3), c(5,4))
    expect_equal(pwDist3[1,1], 0.1551)
    expect_equal(pwDist3[1,2], 0.07132)
    expect_equal(pwDist3[2,1], 0.)
    expect_equal(pwDist3[2,2], 0.08556)
    expect_equal(pwDist3[3,1], 0.)
    expect_equal(pwDist3[3,2], 0.)
    
    # with verbose
    msg <- capture_messages(pwDist4 <- pairwiseEMDDist(
        x = exprList,
        rowRange = c(1,5),
        colRange = c(2,5),
        channels = c("FSC-A", "SSC-A"),
        binSize = 0.05,
        minRange = -10,
        maxRange = 10,
        verbose = TRUE))
    
    pwDist4 <- as.matrix(pwDist4)
    
    expect_equal(length(msg), 22)
    expect_equal(msg[5], "Calculating histogram for file 2...\n")
    expect_equal(msg[10], "Loading file 5...\n")
    expect_equal(msg[21], "i = 3; j = 5; sum(dist) = 0.01813\n")
    
    expect_equal(dim(pwDist4), c(5,4))
    expect_equal(pwDist4[1,1], 0.1551)
    expect_equal(pwDist4[1,2], 0.07132)
    expect_equal(pwDist4[2,1], 0.)
    expect_equal(pwDist4[2,2], 0.08556)
    expect_equal(pwDist4[3,1], 0.)
    expect_equal(pwDist4[3,2], 0.)
})

test_that("pairwiseEMDDist works with fs and BiocParallel", {
    
    logDir <- file.path(outputDir, "BiocParallel", "log")
    
    suppressWarnings(dir.create(logDir, recursive = TRUE))
    
    bp <- BiocParallel::SnowParam(workers = 2, log = TRUE, 
                                  logdir = logDir,
                                  progressbar = TRUE)
    
    pwDist <- suppressWarnings(pairwiseEMDDist(
        x = OMIP021Trans,
        channels = c("FSC-A", "SSC-A"),
        binSize = 0.05,
        minRange = -10,
        maxRange = 10,
        BPPARAM = bp))
    
    pwDist <- as.matrix(pwDist)

    expect_equal(dim(pwDist), c(2,2))
    expect_equal(pwDist[1,1], 0.)
    expect_equal(pwDist[1,2], 0.1551)
    expect_equal(pwDist[2,1], 0.1551)
    expect_equal(pwDist[2,2], 0.)
    
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
    
    fsNames <- c("Donor1", "Donor2", paste0("Agg",1:3))
    names(ffList) <- fsNames
    
    fsAll <- as(ffList,"flowSet")
    
    pwDist2 <- suppressWarnings(pairwiseEMDDist(
            x = fsAll,
            rowRange = c(1,3),
            colRange = c(4,5),
            channels = c("FSC-A", "SSC-A"),
            binSize = 0.05,
            minRange = -10,
            maxRange = 10,
            BPPARAM = bp))
    
    pwDist2 <- as.matrix(pwDist2)

    expect_equal(dim(pwDist2), c(3,2))
    expect_equal(pwDist2[1,1], 0.07241)
    expect_equal(pwDist2[1,2], 0.08347)
    expect_equal(pwDist2[2,1], 0.08501)
    expect_equal(pwDist2[2,2], 0.07451)
    expect_equal(pwDist2[3,1], 0.01293)
    expect_equal(pwDist2[3,2], 0.01813)

    pwDist3 <- suppressWarnings(pairwiseEMDDist(
        x = fsAll,
        rowRange = c(1,2),
        colRange = c(4,5),
        channels = c("FSC-A", "SSC-A"),
        binSize = 0.05,
        minRange = -10,
        maxRange = 10,
        BPPARAM = bp))
    
    pwDist3 <- as.matrix(pwDist3)

    expect_equal(dim(pwDist3), c(2,2))
    expect_equal(pwDist3[1,1], 0.07241)
    expect_equal(pwDist3[1,2], 0.08347)
    expect_equal(pwDist3[2,1], 0.08501)
    expect_equal(pwDist3[2,2], 0.07451)
})

test_that("pairwiseEMDDist works with expr matrix and BiocParallel", {
    
    expr1 <- flowCore::exprs(OMIP021Trans[[1]])[, 
        CytoPipeline::areSignalCols(OMIP021Trans[[1]])]
    expr2 <- flowCore::exprs(OMIP021Trans[[2]])[, 
        CytoPipeline::areSignalCols(OMIP021Trans[[2]])]
    
    logDir <- file.path(outputDir, "BiocParallel", "log")
    
    suppressWarnings(dir.create(logDir, recursive = TRUE))
    
    bp <- BiocParallel::SnowParam(workers = 2, log = TRUE, 
                                  logdir = logDir,
                                  progressbar = TRUE)
    
    pwDist <- suppressWarnings(pairwiseEMDDist(
        x = list(expr1, expr2),
        channels = c("FSC-A", "SSC-A"),
        binSize = 0.05,
        minRange = -10,
        maxRange = 10,
        BPPARAM = bp))
    
    pwDist <- as.matrix(pwDist)
    
    expect_equal(dim(pwDist), c(2,2))
    expect_equal(pwDist[1,1], 0.)
    expect_equal(pwDist[1,2], 0.1551)
    expect_equal(pwDist[2,1], 0.1551)
    expect_equal(pwDist[2,2], 0.)
    
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
    
    fsNames <- c("Donor1", "Donor2", paste0("Agg",1:3))
    names(ffList) <- fsNames
    
    exprList <- lapply(ffList, FUN = function(ff) {
        expr <- flowCore::exprs(ff)[, CytoPipeline::areSignalCols(ff)]
    })
    
    pwDist2 <- suppressWarnings(pairwiseEMDDist(
        x = exprList,
        rowRange = c(1,3),
        colRange = c(4,5),
        channels = c("FSC-A", "SSC-A"),
        binSize = 0.05,
        minRange = -10,
        maxRange = 10,
        BPPARAM = bp))
    
    pwDist2 <- as.matrix(pwDist2)
    
    expect_equal(dim(pwDist2), c(3,2))
    expect_equal(pwDist2[1,1], 0.07241)
    expect_equal(pwDist2[1,2], 0.08347)
    expect_equal(pwDist2[2,1], 0.08501)
    expect_equal(pwDist2[2,2], 0.07451)
    expect_equal(pwDist2[3,1], 0.01293)
    expect_equal(pwDist2[3,2], 0.01813)
    
    pwDist3 <- suppressWarnings(pairwiseEMDDist(
        x = exprList,
        rowRange = c(1,2),
        colRange = c(4,5),
        channels = c("FSC-A", "SSC-A"),
        binSize = 0.05,
        minRange = -10,
        maxRange = 10,
        BPPARAM = bp))
    
    pwDist3 <- as.matrix(pwDist3)
    
    expect_equal(dim(pwDist3), c(2,2))
    expect_equal(pwDist3[1,1], 0.07241)
    expect_equal(pwDist3[1,2], 0.08347)
    expect_equal(pwDist3[2,1], 0.08501)
    expect_equal(pwDist3[2,2], 0.07451)
})

test_that("pairwiseEMDDist dynamic memory loading simulation", {
    datasetPath <- system.file("extdata",
                               package = "CytoPipeline")
    
    files <- list.files(datasetPath, pattern = "Donor", recursive = TRUE)
    if (length(files) != 2) stop("unexpected nb of files!")
    
    filePaths <- file.path(datasetPath, files)
    
    nSamples <- 10
    verbose <- FALSE
    pwDist <- pairwiseEMDDist(
        x = nSamples,
        loadExprMatrixFUN = simulMemoryLoad,
        loadExprMatrixFUNArgs = list(
            filePaths = filePaths,
            nSamples = nSamples,
            verbose = FALSE
        ),
        channels = c("FSC-A", "SSC-A"),
        verbose = FALSE)
    
    pwDist <- as.matrix(pwDist)
    
    expect_equal(dim(pwDist), c(10,10))
    expect_equal(pwDist[1,2], 0.07070)
    expect_equal(pwDist[1,7], 0.14340)
    expect_equal(pwDist[4,8], 0.16625)
    expect_equal(pwDist[6,10], 0.05840)
    
    # same with BiocParallel
    
    logDir <- file.path(outputDir, "BiocParallel", "log")
    
    suppressWarnings(dir.create(logDir, recursive = TRUE))
    bp <- BiocParallel::SnowParam(log = TRUE, 
                                  logdir = logDir,
                                  progressbar = TRUE,
                                  RNGseed = 0)
    verbose <- FALSE
    pwDist2 <- suppressWarnings(pairwiseEMDDist(
        x = nSamples,
        loadExprMatrixFUN = simulMemoryLoad,
        loadExprMatrixFUNArgs = list(
            filePaths = filePaths,
            nSamples = nSamples,
            verbose = verbose
        ),
        channels = c("FSC-A", "SSC-A"),
        verbose = verbose,
        BPPARAM = bp,
        BPOPTIONS = BiocParallel::bpoptions(
            packages = c("flowCore", "CytoPipeline"))))
    
    pwDist2 <- as.matrix(pwDist2)
    
    # Note it is normal that the computed distances are different 
    # from the ones obtained when not using BiocParallel.
    # As described in section 2.4 of technical note:
    # https://bioconductor.org/packages/release/bioc/vignettes/
    # BiocParallel/inst/doc/Random_Numbers.html ,
    # it is not possible to reconcile results from lapply() with
    # results from bplapply().
    # However, using 'RNGseed' argument in BiocParallel::bp() insures
    # results obtained with bplapply() are still reproducible from
    # one run to another
    # 
    expect_equal(dim(pwDist2), c(10,10))
    expect_equal(pwDist2[1,2], 0.07070)
    expect_equal(pwDist2[1,7], 0.14340)
    expect_equal(pwDist2[4,8], 0.16625)
    expect_equal(pwDist2[6,10], 0.05840)
})

test_that("channelSummaryStats works with fs", {

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

    channelsOrMarkers <- c("FSC-A", "SSC-A", "BV785 - CD3")

    ret <- channelSummaryStats(
        fsAll,
        channels = channelsOrMarkers,
        statFUNs = stats::median,
        verbose = FALSE)

    expect_equal(unname(rownames(ret)), flowCore::sampleNames(fsAll))
    expect_equal(unname(colnames(ret)), channelsOrMarkers)

    expect_equal(unname(ret[2,1]), 1.941572819633120)
    expect_equal(unname(ret[3,2]), 1.131832378296435)
    expect_equal(unname(ret[4,3]), 1.638696159955064)

    ret <- channelSummaryStats(
        fsAll,
        channels = channelsOrMarkers,
        statFUNs = list("mean" = mean, "std.dev" = stats::sd),
        verbose = FALSE)

    expect_equal(names(ret), c("mean", "std.dev"))
    expect_equal(unname(rownames(ret[[1]])), flowCore::sampleNames(fsAll))
    expect_equal(unname(colnames(ret[[1]])), channelsOrMarkers)
    expect_equal(unname(ret[[1]][2,1]), 1.961771149533659)
    expect_equal(unname(ret[[1]][3,2]), 1.393034241892733)
    expect_equal(unname(ret[[1]][4,3]), 1.907561024702794)

    expect_equal(unname(rownames(ret[[2]])), flowCore::sampleNames(fsAll))
    expect_equal(unname(colnames(ret[[2]])), channelsOrMarkers)
    expect_equal(unname(ret[[2]][2,1]), 0.5942291770870205)
    expect_equal(unname(ret[[2]][3,2]), 0.7352746696905406)
    expect_equal(unname(ret[[2]][4,3]), 0.8420740225208073)

    # test with only one flow frame
    ret <- channelSummaryStats(
        fsAll[1],
        channels = channelsOrMarkers,
        statFUNs = list("mean" = mean, "std.dev" = stats::sd),
        verbose = FALSE)

    expect_equal(names(ret), c("mean", "std.dev"))
    expect_equal(unname(rownames(ret[[1]])), "Donor1")
    expect_equal(unname(rownames(ret[[2]])), "Donor1")
    expect_equal(unname(colnames(ret[[1]])), channelsOrMarkers)
    expect_equal(unname(colnames(ret[[2]])), channelsOrMarkers)
    expect_equal(unname(ret[[1]][1]), 1.900298)
    expect_equal(unname(ret[[1]][2]), 1.39186533)
    expect_equal(unname(ret[[1]][3]), 1.8544648)
    expect_equal(unname(ret[[2]][1]), 0.73461306)
    expect_equal(unname(ret[[2]][2]), 0.74513873)
    expect_equal(unname(ret[[2]][3]), 0.85498796)

    # one flow frame, one stat
    ret <- channelSummaryStats(
        fsAll[1],
        channels = channelsOrMarkers,
        statFUNs = list("mean" = mean),
        verbose = FALSE)

    expect_equal(unname(rownames(ret[[1]])), "Donor1")
    expect_equal(unname(colnames(ret[[1]])), channelsOrMarkers)
    expect_equal(unname(ret[[1]][1,1]), 1.900298)
    expect_equal(unname(ret[[1]][1,2]), 1.39186533)
    expect_equal(unname(ret[[1]][1,3]), 1.8544648)

    # one flow frame, one stat (bis with direct impact of stat FUN)
    ret <- channelSummaryStats(
        fsAll[1],
        channels = channelsOrMarkers,
        statFUNs = mean,
        verbose = FALSE)

    expect_equal(unname(rownames(ret)), "Donor1")
    expect_equal(unname(colnames(ret)), channelsOrMarkers)
    expect_equal(unname(ret[1,1]), 1.900298)
    expect_equal(unname(ret[1,2]), 1.39186533)
    expect_equal(unname(ret[1,3]), 1.8544648)

    # case where no channels is provided
    ret <- channelSummaryStats(
        fsAll,
        statFUNs = list("mean" = mean, "std.dev" = stats::sd),
        verbose = FALSE
    )

    allSignalChannels <-
        flowCore::colnames(fsAll)[CytoPipeline::areSignalCols(fsAll)]
    nSignalCh <- length(allSignalChannels)

    allSignalChannelNames <- allSignalChannels
    for (i in seq_along(allSignalChannels)) {
        channelMarker <-
            flowCore::getChannelMarker(fsAll[[1]], allSignalChannels[i])$desc
        if (!is.null(channelMarker) && !is.na(channelMarker)){
            allSignalChannelNames[i] <- channelMarker
        }
    }

    expect_equal(unname(rownames(ret[[1]])),
                 c("Donor1", "Donor2", "Agg1", "Agg2", "Agg3"))
    expect_equal(unname(colnames(ret[[1]])), allSignalChannelNames)

})

test_that("channelSummaryStats works with expr matrix", {
    
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
    
    exprList <- lapply(
        ffList,
        FUN = function(ff) {
            exprMat <- 
                flowCore::exprs(ff)[,CytoPipeline::areSignalCols(ff)]
            exprMat
        }
    )
    
    names(exprList) <- names(ffList)
    
    channels <- c("FSC-A", "SSC-A", "670/30Violet-A")
    
    ret <- channelSummaryStats(
        exprList,
        channels = channels,
        statFUNs = stats::median,
        verbose = FALSE)
    
    expect_equal(unname(rownames(ret)), names(ffList))
    expect_equal(unname(colnames(ret)), channels)
    
    expect_equal(unname(ret[2,1]), 1.941572819633120)
    expect_equal(unname(ret[3,2]), 1.131832378296435)
    expect_equal(unname(ret[4,3]), 1.638696159955064)
    
    ret <- channelSummaryStats(
        exprList,
        channels = channels,
        statFUNs = list("mean" = mean, "std.dev" = stats::sd),
        verbose = FALSE)
    
    expect_equal(names(ret), c("mean", "std.dev"))
    expect_equal(unname(rownames(ret[[1]])), names(ffList))
    expect_equal(unname(colnames(ret[[1]])), channels)
    expect_equal(unname(ret[[1]][2,1]), 1.961771149533659)
    expect_equal(unname(ret[[1]][3,2]), 1.393034241892733)
    expect_equal(unname(ret[[1]][4,3]), 1.907561024702794)
    
    expect_equal(unname(rownames(ret[[2]])), names(ffList))
    expect_equal(unname(colnames(ret[[2]])), channels)
    expect_equal(unname(ret[[2]][2,1]), 0.5942291770870205)
    expect_equal(unname(ret[[2]][3,2]), 0.7352746696905406)
    expect_equal(unname(ret[[2]][4,3]), 0.8420740225208073)
    
    # test with only one flow frame
    ret <- channelSummaryStats(
        exprList[1],
        channels = channels,
        statFUNs = list("mean" = mean, "std.dev" = stats::sd),
        verbose = FALSE)
    
    expect_equal(names(ret), c("mean", "std.dev"))
    expect_equal(unname(rownames(ret[[1]])), "Donor1")
    expect_equal(unname(rownames(ret[[2]])), "Donor1")
    expect_equal(unname(colnames(ret[[1]])), channels)
    expect_equal(unname(colnames(ret[[2]])), channels)
    expect_equal(unname(ret[[1]][1]), 1.900298)
    expect_equal(unname(ret[[1]][2]), 1.39186533)
    expect_equal(unname(ret[[1]][3]), 1.8544648)
    expect_equal(unname(ret[[2]][1]), 0.73461306)
    expect_equal(unname(ret[[2]][2]), 0.74513873)
    expect_equal(unname(ret[[2]][3]), 0.85498796)
    
    # one flow frame, one stat
    ret <- channelSummaryStats(
        exprList[1],
        channels = channels,
        statFUNs = list("mean" = mean),
        verbose = FALSE)
    
    expect_equal(unname(rownames(ret[[1]])), "Donor1")
    expect_equal(unname(colnames(ret[[1]])), channels)
    expect_equal(unname(ret[[1]][1,1]), 1.900298)
    expect_equal(unname(ret[[1]][1,2]), 1.39186533)
    expect_equal(unname(ret[[1]][1,3]), 1.8544648)
    
    # one flow frame, one stat (bis with direct impact of stat FUN)
    ret <- channelSummaryStats(
        exprList[1],
        channels = channels,
        statFUNs = mean,
        verbose = FALSE)
    
    expect_equal(unname(rownames(ret)), "Donor1")
    expect_equal(unname(colnames(ret)), channels)
    expect_equal(unname(ret[1,1]), 1.900298)
    expect_equal(unname(ret[1,2]), 1.39186533)
    expect_equal(unname(ret[1,3]), 1.8544648)
    
    # case where no channels is provided
    ret <- channelSummaryStats(
        exprList,
        statFUNs = list("mean" = mean, "std.dev" = stats::sd),
        verbose = FALSE
    )
    
    allSignalChannels <- colnames(exprList[[1]])
    nSignalCh <- length(allSignalChannels)
    
    
    expect_equal(unname(rownames(ret[[1]])),
                 c("Donor1", "Donor2", "Agg1", "Agg2", "Agg3"))
    expect_equal(unname(colnames(ret[[1]])), allSignalChannels)
    
})

test_that("channelSummaryStats dynamic memory loading simulation", {
    datasetPath <- system.file("extdata",
                               package = "CytoPipeline")

    files <- list.files(datasetPath, pattern = "Donor", recursive = TRUE)
    if (length(files) != 2) stop("unexpected nb of files!")

    filePaths <- file.path(datasetPath, files)

    channelsOrMarkers <- c("FSC-A", "SSC-A", "BV785 - CD3")

    nSamples <- 10
    verbose <- FALSE
    ret <- channelSummaryStats(
        x = nSamples,
        loadExprMatrixFUN = simulMemoryLoad,
        loadExprMatrixFUNArgs = list(
            filePaths = filePaths,
            nSamples = nSamples,
            channels = channelsOrMarkers,
            verbose = FALSE
        ),
        statFUNs = list("mean" = mean, "std.dev" = stats::sd),
        verbose = verbose)

    expect_equal(names(ret), c("mean", "std.dev"))
    expect_equal(dim(ret[[1]]), c(10, 3))
    expect_equal(unname(colnames(ret[[1]])), channelsOrMarkers)
    expect_equal(unname(colnames(ret[[2]])), channelsOrMarkers)
    expect_equal(unname(ret[[1]][1,1]), 1.90567984)
    expect_equal(unname(ret[[1]][2,2]), 1.36523328)
    expect_equal(unname(ret[[1]][7,3]), 1.99832484)
    expect_equal(unname(ret[[2]][1,1]), 0.71358848)
    expect_equal(unname(ret[[2]][2,2]), 0.74298778)
    expect_equal(unname(ret[[2]][7,3]), 0.82220517)

    # same with BiocParallel

    logDir <- file.path(outputDir, "BiocParallel", "log")

    suppressWarnings(dir.create(logDir, recursive = TRUE))
    bp <- BiocParallel::SnowParam(log = TRUE,
                                  logdir = logDir,
                                  progressbar = TRUE,
                                  RNGseed = 0)
    ret2 <- suppressWarnings(CytoMDS::channelSummaryStats(
        x = nSamples,
        loadExprMatrixFUN = simulMemoryLoad,
        loadExprMatrixFUNArgs = list(
            filePaths = filePaths,
            nSamples = nSamples,
            channels = channelsOrMarkers,
            verbose = FALSE
        ),
        statFUNs = list("mean" = mean, "std.dev" = stats::sd),
        verbose = verbose,
        BPPARAM = bp,
        BPOPTIONS = BiocParallel::bpoptions(
            packages = c("flowCore", "CytoPipeline"))))

    # Note it is normal that the computed stats are different
    # from the ones obtained when not using BiocParallel.
    # As described in section 2.4 of technical note:
    # https://bioconductor.org/packages/release/bioc/vignettes/
    # BiocParallel/inst/doc/Random_Numbers.html ,
    # it is not possible to reconcile results from lapply() with
    # results from bplapply().
    # However, using 'RNGseed' argument in BiocParallel::bp() insures
    # results obtained with bplapply() are still reproducible from
    # one run to another
    #
    expect_equal(names(ret2), c("mean", "std.dev"))
    expect_equal(dim(ret2[[1]]), c(10, 3))
    expect_equal(unname(colnames(ret2[[1]])), channelsOrMarkers)
    expect_equal(unname(colnames(ret2[[2]])), channelsOrMarkers)
    expect_equal(unname(ret2[[1]][1,1]), 1.90567984)
    expect_equal(unname(ret2[[1]][2,2]), 1.36523328)
    expect_equal(unname(ret2[[1]][7,3]), 1.99832484)
    expect_equal(unname(ret2[[2]][1,1]), 0.71358848)
    expect_equal(unname(ret2[[2]][2,2]), 0.74298778)
    expect_equal(unname(ret2[[2]][7,3]), 0.82220517)
})


test_that("computeMetricMDS works", {
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
    
    pwDist <- pairwiseEMDDist(fsAll, 
                              channels = c("FSC-A", "SSC-A"),
                              verbose = FALSE)
    
    mdsObj <- computeMetricMDS(pwDist, nDim = 2, seed = 0)
    
    expect_equal(stress(mdsObj), 0.0203635387)
    tgtspp <- c(
        11.111867, 13.086625, 4.448069, 37.334795, 34.018645)
    names(tgtspp) <- 1:5
    expect_equal(spp(mdsObj), tgtspp)
    
    # do the same but with the full distance object, 
    # and select channels in computeMetricMDS
    pwDistFull <- pairwiseEMDDist(fsAll, 
                                  channels = NULL,
                                  verbose = TRUE)
    
    mdsObjFull <- computeMetricMDS(pwDistFull, nDim = 2, seed = 0,
                                   whichChannels = c("FSC-A", "SSC-A"))
    expect_equal(stress(mdsObjFull), 0.0203635387)
    expect_equal(spp(mdsObjFull), tgtspp)
    
    
    # with no user provided nDim, but (implicit) target pseudo rsquare = 0.95
    mdsObj2 <- computeMetricMDS(pwDist, seed = 0)
    
    expect_equal(nDim(mdsObj2), 2)
    expect_equal(RSqVec(mdsObj2)[2], 0.99843722)
    
    # with no user provided nDim, but explicit target pseudo rsquare = 0.98
    mdsObj3 <- computeMetricMDS(pwDist, seed = 0, targetPseudoRSq = 0.98)
    
    expect_equal(nDim(mdsObj3), 2)
    expect_equal(RSqVec(mdsObj3)[2], 0.99843722)
    
    # with no user provided nDim, but explicit target pseudo rsquare = 0.999
    mdsObj4 <- computeMetricMDS(pwDist, seed = 0, targetPseudoRSq = 0.999)
    
    expect_equal(nDim(mdsObj4), 3)
    expect_equal(RSqVec(mdsObj4)[3], 0.99988906)
    
})

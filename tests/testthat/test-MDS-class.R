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

nHD <- 10
nLD <- 2
nPoints <- 20 
 
# generate uniformly distributed points in 10 dimensions
set.seed(0)
points <- matrix(
     data = runif(n = nPoints * nHD),
     nrow = nPoints)
     
# calculate euclidian distances     
pwDist  <- as.matrix(dist(points))

# compute Metric MDS object by reaching a target pseudo RSquare
mdsObj <- computeMetricMDS(pwDist, targetPseudoRSq = 0.95, seed = 0)

test_that("basic MDS class works", {
    ret <- validObject(mdsObj)
    expect_true(ret)
    
    nDim <- nDim(mdsObj)
    expect_equal(nDim, 7)
    
    pwDist <- pwDist(mdsObj)
    expect_equal(dim(pwDist), c(20,20))
    expect_equal(as.vector(pwDist)[6], 1.2781415)
    
    proj <- unname(projections(mdsObj))
    expect_equal(dim(proj), c(20, nDim))
    expect_equal(proj[13,3], 0.230974185)
    
    projDist <- projDist(mdsObj)
    expect_equal(dim(projDist), c(20,20))
    expect_equal(as.vector(projDist)[7], 1.0551025)
    
    eigen <- eigenVals(mdsObj)
    expect_equal(length(eigen), nDim)
    expect_equal(eigen[2], 2.1465394)
    
    pctvar <- pctvar(mdsObj)
    expect_equal(length(pctvar), nDim)
    expect_equal(pctvar[1], 0.27773307)
    
    RSq <- RSq(mdsObj)
    expect_equal(RSq, 0.98236885)
    
    RSqVec <- RSqVec(mdsObj)
    expect_equal(length(RSqVec), nDim)
    expect_equal(RSqVec[nDim], RSq)
    
    GoF <- GoF(mdsObj)
    expect_equal(length(GoF), nDim)
    expect_equal(GoF[nDim], 0.99939197)
    
    spp <- unname(spp(mdsObj))
    expect_equal(length(spp), nPoints)
    expect_equal(spp[17], 1.21752488)
    
    stress <- stress(mdsObj)
    expect_equal(stress, 0.024658166)
    smacofRes <- smacofRes(mdsObj)
    expect_equal(class(smacofRes), c("smacofB", "smacof"))
})

test_that("MDS object export and reimport works", {
    outputFile <- base::tempfile()
    saveRDS(object = mdsObj, file = outputFile)
    mdsObj2 <- readRDS(file = outputFile)
    ret <- validObject(mdsObj2)
    expect_true(ret)
})
    
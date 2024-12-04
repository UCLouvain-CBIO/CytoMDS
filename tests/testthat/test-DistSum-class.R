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

set.seed(20241125)
nDistr <- 5
nFeat <- 7
M <- matrix(data = rnorm(nDistr * nFeat), ncol = nFeat)
rownames(M) <- paste0("distr", 1:nDistr)
colnames(M) <- paste0("feat", 1:nFeat)

DList <- lapply(colnames(M),
                FUN = function(colName) {
                    D <- as.matrix(dist(
                        M[, colName, drop = FALSE]))
                    D
                })

D <- Reduce(x = DList, f = function(A, B) A + B)

names(DList) <- colnames(M)

test_that("DistSum built with matrix object works", {
    distObj <- DistSum(D)
    #expect_no_error(show(distObj))
    expect_equal(dim(distObj), c(nDistr, nDistr))
    expect_equal(ncol(distObj), nDistr)
    expect_equal(nrow(distObj), nDistr)
    expect_equal(nFeatures(distObj), 1)
    expect_equal(dimnames(distObj)[[1]], rownames(M))
    expect_equal(dimnames(distObj)[[2]], rownames(M))
    expect_equal(rownames(distObj), rownames(M))
    expect_equal(colnames(distObj), rownames(M))
    expect_null(featureNames(distObj))
    dd <- as.matrix(distObj)
    expect_equal(dd, D)

    dd2 <- as.matrix(distObj, whichFeatures = 1)
    expect_equal(dd2, dd)

    expect_error(featureNames(distObj) <- c("ft1", "ft2"),
                 regexp = "new feature names should have same length")
    featureNames(distObj) <- "FT1"
    expect_equal(featureNames(distObj), "FT1")

    expect_error(colnames(distObj) <- c("D1", "D2"),
                 regexp = "not equal to array extent")
    colnames(distObj) <- paste0("D", 1:nDistr)
    expect_equal(colnames(distObj), paste0("D", 1:nDistr))
    
    expect_error(rownames(distObj) <- c("D1", "D2"),
                 regexp = "not equal to array extent")
    rownames(distObj) <- paste0("D", 1:nDistr)
    expect_equal(rownames(distObj), paste0("D", 1:nDistr))
    
    expect_error(dimnames(distObj) <- paste0("D", 1:nDistr),
                 regexp = "must be a list")
    expect_error(dimnames(distObj) <- list(c("D1", "D2")),
                 regexp = "not equal to array extent")
    dimnames(distObj) <- list(paste0("D", 1:nDistr),
                              paste0("D", 1:nDistr))
    expect_equal(dimnames(distObj)[[1]], paste0("D", 1:nDistr))
    expect_equal(dimnames(distObj)[[2]], paste0("D", 1:nDistr))
})

test_that("DistSum built with list works", {
    distObj <- DistSum(DList)
    
    #expect_no_error(show(distObj))
    expect_equal(dim(distObj), c(nDistr, nDistr))
    expect_equal(ncol(distObj), nDistr)
    expect_equal(nrow(distObj), nDistr)
    expect_equal(nFeatures(distObj), nFeat)
    expect_equal(dimnames(distObj)[[1]], rownames(M))
    expect_equal(dimnames(distObj)[[2]], rownames(M))
    expect_equal(rownames(distObj), rownames(M))
    expect_equal(colnames(distObj), rownames(M))
    expect_equal(featureNames(distObj), names(DList))
    
    expect_error(featureNames(distObj) <- c("ft1", "ft2"),
                 regexp = "new feature names should have same length")
    featureNames(distObj) <- paste0("FT", seq(nFeat))
    expect_equal(featureNames(distObj), paste0("FT", seq(nFeat)))
    
    expect_error(colnames(distObj) <- c("D1", "D2"),
                 regexp = "not equal to array extent")
    colnames(distObj) <- paste0("D", 1:nDistr)
    expect_equal(colnames(distObj), paste0("D", 1:nDistr))
    
    expect_error(rownames(distObj) <- c("D1", "D2"),
                 regexp = "not equal to array extent")
    rownames(distObj) <- paste0("D", 1:nDistr)
    expect_equal(rownames(distObj), paste0("D", 1:nDistr))
    
    expect_error(dimnames(distObj) <- paste0("D", 1:nDistr),
                 regexp = "must be a list")
    expect_error(dimnames(distObj) <- list(c("D1", "D2")),
                 regexp = "not equal to array extent")
    dimnames(distObj) <- list(paste0("D", 1:nDistr), 
                              paste0("D", 1:nDistr))
    expect_equal(dimnames(distObj)[[1]], paste0("D", 1:nDistr))
    expect_equal(dimnames(distObj)[[2]], paste0("D", 1:nDistr))
    
    # reset everything
    distObj <- DistSum(DList)
    
    dd <- as.matrix(distObj)
    expect_equal(dd, D)
    expect_equal(as.matrix(dd)[1,2], 9.58176716)
    
    dd1 <- as.matrix(distObj, whichFeatures = 1)
    expect_equal(as.matrix(dd1)[4,5],0.02922822)
    
    dd1bis <- as.matrix(distObj, whichFeatures = colnames(M)[1])
    expect_equal(dd1bis, dd1)
    
    ddList <- lapply(colnames(M),
                     function(colname){
                         as.matrix(distObj, whichFeatures = colname)
                     })
    
    ddSum <- Reduce(x = ddList,
                    f = function(x, y) x+y)
    
    expect_equal(dd, ddSum)
    
    
    ddPart <- as.matrix(distObj, whichFeatures = colnames(M)[1:2])
    dd2 <- as.matrix(distObj, whichFeatures = 2)
    
    expect_equal(ddPart, dd1+dd2)
    

})

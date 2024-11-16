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

#' @title MDS class
#' 
#' @rdname MDS
#'
#' @description
#'
#' Class representing Multi Dimensional Scaling (MDS) projection.  
#' 
#' @slot nDim `numeric`, nb of dimensions of the projection
#'
#' @slot pwDist An object of class `dist` storing 
#' the triangular relevant part of the symmetric, zero diagonal 
#' pairwise distance matrix (nPoints * nPoints), BEFORE projection.  
#'
#' @slot proj The projection matrix, resulting from MDS  
#' 
#' @slot projDist An object of class `dist` storing 
#' the triangular relevant part of the symmetric, zero diagonal 
#' pairwise distance matrix (nPoints * nPoints), AFTER projection. 
#'
#' @slot eigen `numeric`, vector of `nDim` length, containing the eigen 
#' values of the PCA that is applied after the Smacof algorithm. 
#'
#' @slot pctvar `numeric`, vector of `nDim` length, containing the percentage 
#' of explained variance per axis.  
#'
#' @slot RSq `numeric`, vector of pseudo R square indicators, 
#' as a function of number of dimensions. 
#' `RSq[nDim]` is the global pseudo R square, as displayed on plots. 
#' 
#' @slot GoF `numeric`, vector of goodness of fit indicators,
#' as a function of number of dimensions.
#' `GoF[nDim]` is the global goodness of fit.   
#' 
#' Note pseudo R square and goodness of fit indicators are essentially the 
#' same indicator, only the definition of total sum of squares differ:
#' - for pseudo RSq: TSS is calculated using the mean pairwise distance 
#' as minimum
#' - for goodness of fit: TSS is calculated using 0 as minimum
#' 
#' @slot smacofRes an object of class 'smacofB' containing the algorithmic 
#' optimization results, for example stress and stress per point, 
#' as returned by `smacof::smacofSym()` method.
#' 
#' @exportClass MDS
#' 
#' @return nothing
#' 
#' @examples
#' 
#' 
#' nHD <- 10
#' nLD <- 2
#' nPoints <- 20 
#' 
#' # generate uniformly distributed points in 10 dimensions
#' points <- matrix(
#'     data = runif(n = nPoints * nHD),
#'     nrow = nPoints)
#'     
#' # calculate euclidian distances     
#' pwDist  <- dist(points)
#' 
#' # compute Metric MDS object by reaching a target pseudo RSquare
#' mdsObj <- computeMetricMDS(pwDist, targetPseudoRSq = 0.95)
#' 
#' show(mdsObj)
#' 
#'
setClass("MDS",
         slots = c(
             nDim = "numeric",
             pwDist = "dist",
             proj = "matrix",
             projDist = "dist",
             eigen = "numeric",
             pctvar = "numeric",
             RSq = "numeric",
             GoF = "numeric",
             smacofRes = "ANY"
         ),
         prototype = list(
             nDim = 0,
             pwDist = dist(numeric()),
             proj = matrix(),
             projDist = dist(numeric()),
             eigen = numeric(),
             pctvar = numeric(),
             RSq = numeric(),
             GoF = numeric(),
             smacofRes = NULL
         )
)

setValidity("MDS", function(object) {
    if (object@nDim > 0){
        nPoints <- nPoints(object)
        if (nrow(object@proj) != nPoints){
            return("proj matrix: row number inconsistency")
        }
        if (ncol(object@proj) != object@nDim){
            return("proj matrix: col number inconsistency")
        }
        if (!all.equal(dim(object@projDist), dim(object@pwDist))) {
            return("proj distances: dimension inconsistency")
        }
        if (length(object@eigen) != object@nDim) {
            return("eigen vector: length inconsistency")
        }
        if (length(object@pctvar) != object@nDim) {
            return("pctvar vector: length inconsistency")
        }
        if (length(object@RSq) != object@nDim) {
            return("RSq vector: length inconsistency")
        }
        if (length(object@GoF) != object@nDim) {
            return("GoF vector: length inconsistency")
        }
        if (is.null(object@smacofRes)) {
            return("Un-populated smacofRes")
        }
        if (!inherits(object@smacofRes, "smacofB")) {
            return("smacofRes should inherit from smacofB class")
        }
    }
    
    return(TRUE)

})

#' @rdname MDS
#' @param object a `MDS` object
#'
setMethod(
    "show", "MDS",
    function(object) {
        cat("MDS object containing MDS projection (using Smacof algorithm) ",
            "data:\n")
        cat("Nb of dimensions: ", object@nDim, "\n")
        if (object@nDim > 0) {
            cat("Nb of points: ", nrow(object@proj), "\n")
            cat("Stress: ", round(object@smacofRes$stress, 6), "\n")
            cat("Pseudo RSquare: ", round(object@RSq[object@nDim], 6), "\n")
            cat("Goodness of fit: ", round(object@GoF[object@nDim], 6), "\n")
        }
    }
)

#' @rdname MDS
#' @param x a `MDS` object
#' @export
nDim <- function(x) {
    stopifnot(inherits(x, "MDS"))
    return(x@nDim)
}

#' @rdname MDS
#' @param x a `MDS` object
#' @export
nPoints <- function(x) {
    stopifnot(inherits(x, "MDS"))
    return(attr(x@pwDist, "Size"))
}

#' @rdname MDS
#' @param x a `MDS` object
#' @export
pwDist <- function(x) {
    stopifnot(inherits(x, "MDS"))
    return(x@pwDist)
}

#' @rdname MDS
#' @param x a `MDS` object
#' @export
projections <- function(x) {
    stopifnot(inherits(x, "MDS"))
    return(x@proj)
}

#' @rdname MDS
#' @param x a `MDS` object
#' @export
projDist <- function(x) {
    stopifnot(inherits(x, "MDS"))
    return(x@projDist)
}

#' @rdname MDS
#' @description
#' returns the value of the stress criterion, minimized by the 
#' SMACOF algorithm.
#' @param x a `MDS` object
#' @export
stress <- function(x) {
    stopifnot(inherits(x, "MDS"))
    return(x@smacofRes$stress)
}

#' @rdname MDS
#' @description
#' returns a vector of nPoints dimension, containing the stress 
#' indicator per point. The `stress` minimization criterion can indeed be 
#' allocated per represented point. The more the stress of a particular point, 
#' the less accurate its distances w.r.t. the other points.
#' @param x a `MDS` object
#' @export
spp <- function(x) {
    stopifnot(inherits(x, "MDS"))
    return(x@smacofRes$spp)
}

#' @rdname MDS
#' @param x a `MDS` object
#' @export
eigenVals <- function(x) {
    stopifnot(inherits(x, "MDS"))
    return(x@eigen)
}

#' @rdname MDS
#' @param x a `MDS` object
#' @export
pctvar <- function(x) {
    stopifnot(inherits(x, "MDS"))
    return(x@pctvar)
}

#' @rdname MDS
#' @param x a `MDS` object
#' @export
RSq <- function(x) {
    stopifnot(inherits(x, "MDS"))
    return(x@RSq[x@nDim])
}

#' @rdname MDS
#' @param x a `MDS` object
#' @export
RSqVec <- function(x) {
    stopifnot(inherits(x, "MDS"))
    return(x@RSq)
}

#' @rdname MDS
#' @param x a `MDS` object
#' @export
GoF <- function(x) {
    stopifnot(inherits(x, "MDS"))
    return(x@GoF)
}

#' @rdname MDS
#' @param x a `MDS` object
#' @export
smacofRes <- function(x) {
    stopifnot(inherits(x, "MDS"))
    return(x@smacofRes)
}


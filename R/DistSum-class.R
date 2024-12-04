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

#' @title DistSum class
#' 
#' @rdname DistSum-class
#' 
#' @name DistSum
#' 
#' @aliases DistSum-class
#'
#' @description
#'
#' Class representing pairwise distances between multiple multidimensional 
#' distributions, when the distance is calculated as a sum of marginal 
#' distribution distances.
#' 
#' @slot pwFullDist A `matrix`, storing the full pairwise distance matrix, 
#' calculated as the sum of the distances per dimension (feature). 
#' Note this matrix is not necessarily a square symmetric matrix, 
#' as it could be ocasionally used to store blocks of the distance matrix. 
#' 
#' @slot pwDistPerDim A `list` of `matrix` objects storing the contribution
#' of each dimension (feature) of the multidimensional distributions
#' to the full distance matrix
#'
#' @exportClass DistSum
#' 
#' @return nothing
#' @examples
#' 
#' # create a dummy distance matrix 
#' # to do this we use `nPoints` points 
#' # in an euclidian space of `nFeat` dimensions
#' nPoints <- 5
#' nFeat <- 7
#' M <- matrix(data = rnorm(nPoints * nFeat), ncol = nFeat)
#' rownames(M) <- paste0("point", 1:nPoints)
#' colnames(M) <- paste0("feat", 1:nFeat)
#' 
#' DList <- lapply(colnames(M),
#' FUN = function(colName) {
#'     D <- as.matrix(dist(
#'         M[, colName, drop = FALSE]))
#'     D
#' })
#' 
#' D <- Reduce(x = DList, f = function(A, B) A + B)
#' 
#' names(DList) <- colnames(M)
#' 
#' # Example of creating of a DistSum object based on the full distance matrix
#' distObj1 <- DistSum(D)
#' show(distObj1)
#' 
#' # Example of creation of a DistSum object based on a list of matrices
#' # representing the additive contribution of each feature
#' distObj2 <- DistSum(DList)
#' 
#' show(distObj2)
#'
#' # getting dimensions
#' myDim <- dim(distObj2) # c(nPoints, nPoints)
#' ncols <- ncol(distObj2) # nPoints
#' nrows <- nrow(distObj2) # nPoints
#' nFeats <- nFeatures(distObj2) # nFeat
#' myFeatNames <- featureNames(distObj2) # paste0("feat", 1:nFeat)
#' myRowNames <- rownames(distObj2) # paste0("point", 1:nPoints)
#' myRowNames <- colnames(distObj2) # paste0("point", 1:nPoints)
#' 
#' # get full distance matrix
#' dd <- as.matrix(distObj2)
#' 
#' # get partial distance matrix for feature 1
#' dd1 <- as.matrix(distObj2, whichFeatures = 1)
#' 
#' # same thing, using feature name
#' dd1bis <- as.matrix(distObj2, whichFeatures = "feat1")
#' 
#' # getting partial distance for feature 1 & 2
#' 
#' ddPart <- as.matrix(distObj2, whichFeatures = colnames(M)[1:2])

setClass("DistSum",
         slots = c(
             pwFullDist = "matrix",
             pwDistPerDim = "list"
         ),
         prototype = list(
             pwFullDist = matrix(numeric()),
             pwDistPerDim = list()
         )
)

#' @importFrom methods is
.isValid <- function(object) {
    nc <- ncol(object)
    nr <- nrow(object)
    cnames <- colnames(object)
    rnames <- rownames(object)
    if (!is.null(cnames)) {
        if (length(unique(cnames)) !=
            length(cnames)) {
            return("duplicate column names found")
        }
    }
    if (!is.null(rnames)) {
        if (length(unique(rnames)) !=
            length(rnames)) {
            return("duplicate row names found")
        }
    }
    if (nc * nr > 0) {
        nFeat <- nFeatures(object)
        if (nFeat > 0) {
            featNames <- featureNames(object)
            if (!is.null(featNames)) {
                if (length(unique(featNames)) !=
                    length(featNames)) {
                    return("duplicate feature names found")
                }
                if (length(featNames) != nFeat) {
                    return(paste0("length of featureNames array does not ",
                                  "returned number of features"))
                }
            }
        } else {
            return("no features found in DistSum object")
        }
    }
    
    #check classes of main objects
    if (!is(object@pwFullDist, "matrix")) {
        return("pwFullDist object is not a `matrix` object")
    }
    
    ret <- lapply(
        object@pwDistPerDim,
        FUN = function(item, dims) {
            if (!is(item, "matrix")) {
                return(paste0("pwDistPerDim list contains other objects ",
                              "than `matrix` objects"))
            }
            if (is.null(dims) && !is.null(dim(item))) {
                return(paste0("pwDistPerDim list contains matrices with ",
                              "not null dimensions while pwFullDist ",
                              "has null dimensions"))
            }
            if (!identical(dim(item), dims)) {
                return(paste0("pwDistPerDim list contains matrices with ",
                              "dimensions not equal to pwFullDist dimensions"))
            }
        },
        dims = dim(object@pwDistPerDim))
    
    return(TRUE)
}

setValidity("DistSum", function(object) {
    return(.isValid(object))
})

#' @rdname DistSum-class
#' @param object a `DistSum` object
#'
#' @importMethodsFrom methods show
setMethod(
    "show", "DistSum",
    function(object) {
        cat("`DistSum` object containing pairwise distances", 
            "between distributions\n", 
            "and their decomposition as a sum of feature contributions\n")
        cat("Matrix dimensions: ", dim(object), "\n")
        if (ncol(object) * nrow(object) > 0) {
            cat("Nb of features: ", nFeatures(object), "\n")
            cat("Feature names: ", featureNames(object), "\n")
            cat("Full distance matrix: \n")
            print(.as.matrix(object))
        }
    }
)

setGeneric("DistSum", function(object, ...) {
    standardGeneric("DistSum")
})

#' @rdname DistSum-class
#' @param object a `dist` object containing the full distance matrix
#' @export
setMethod("DistSum",
          "matrix",
          function(object) {
              newObj <- methods::new("DistSum",
                                     pwFullDist = object,
                                     pwDistPerDim = list(object))


              if (isTRUE(methods::validObject(newObj))){
                  return(newObj)
              }
          })

#' @rdname DistSum-class
#' @param object a `list` of which each element is a `dist` object representing
#' the distance matrix for one feature.
#' @export
setMethod(
    "DistSum", "list",
    function(object) {
        pwFUllDist <- Reduce(
            f = function(x, y){
                x + y
            },
            x = object)
        
        newObj <- methods::new("DistSum",
                               pwFullDist = pwFUllDist,
                               pwDistPerDim = object)
        
        if (isTRUE(methods::validObject(newObj))){
            return(newObj)
        }
    }
)

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @export
setMethod(
    "dim", "DistSum",
    function(x) {
        dd <- dim(x@pwFullDist)
        return(dd)
    }
)

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @export
setMethod(
    "dimnames", "DistSum",
    function(x) {
        dn <- dimnames(x@pwFullDist)
        return(dn)
    }
)

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @param value the new dimension names to be assigned
#' @export
setMethod(
    "dimnames<-", c(x = "DistSum", value = "character"),
    function(x, value) {
        dimnames(x@pwFullDist) <- value
        lapply(x@pwDistPerDim,
               FUN = function(mat) {
                   dimnames(mat) <- value
                   mat
               })
        if (isTRUE(methods::validObject(x)))
            return(x)
    }
)

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @param value the new dimension names to be assigned
#' @export
setMethod(
    "dimnames<-", c(x = "DistSum", value = "ANY"),
    function(x, value) {
        dimnames(x@pwFullDist) <- value
        lapply(x@pwDistPerDim,
               FUN = function(mat) {
                  dimnames(mat) <- value
                   mat
               })
        if (isTRUE(methods::validObject(x)))
            return(x)
    }
)

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @export
setMethod(
    "ncol", "DistSum",
    function(x) {
        nc <- ncol(x@pwFullDist)
        return(nc)
    }
)

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @export
setMethod(
    "colnames", "DistSum",
    function(x) {
        cn <- colnames(x@pwFullDist)
        return(cn)
    }
)

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @param value the new column names to be assigned
#' @export
setMethod(
    "colnames<-", "DistSum",
    function(x, value) {
        colnames(x@pwFullDist) <- value
        lapply(x@pwDistPerDim,
               FUN = function(mat) {
                   colnames(mat) <- value
                   mat
               })
        if (isTRUE(methods::validObject(x))) return(x)
    }
)

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @export
setMethod(
    "nrow", "DistSum",
    function(x) {
        nr <- nrow(x@pwFullDist)
        return(nr)
    }
)

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @export
setMethod(
    "rownames", "DistSum",
    function(x) {
        rn <- rownames(x@pwFullDist)
        return(rn)
    }
)

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @param value the new column names to be assigned
#' @export
setMethod(
    "rownames<-", "DistSum",
    function(x, value) {
        rownames(x@pwFullDist) <- value
        lapply(x@pwDistPerDim,
               FUN = function(mat) {
                   rownames(mat) <- value
                   mat
               })
        if (isTRUE(methods::validObject(x))) return(x)
    }
)

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @export
nFeatures <- function(x) {
    stopifnot(inherits(x, "DistSum"))
    return(length(x@pwDistPerDim))
}

#' @rdname DistSum-class
#' @param object a `DistSum` object
#' @importMethodsFrom Biobase featureNames
#' @export
setMethod(
    "featureNames", "DistSum",
    function(object) {
     return(names(object@pwDistPerDim))
    }
)

#' @rdname DistSum-class
#' @param object a `DistSum` object
#' @param value the new feature names to be assigned
#' @importMethodsFrom Biobase featureNames<-
#' @export
setMethod(
    "featureNames<-", "DistSum",
    function(object, value){
        if (!is.character(value)) {
            stop("new feature names should be a character vector")
        }
        if (length(value) != nFeatures(object)) {
            stop("new feature names should have same length as nb of features")
        }
        names(object@pwDistPerDim) <- value
        
        if (isTRUE(methods::validObject(object))) return(object)
    }
)

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @param i the row index
#' @param j the column index
#' @param ... other arguments (not used)
#' @param drop if TRUE, decrease the nb of dimensions when possible
#' @export
setMethod(
    '[', c(x = "DistSum", i = "ANY", j = "ANY", drop = "ANY"),
    function(x, i, j, ..., drop) {
        .as.matrix(x)[i = i, j = j, ..., drop = drop]
    })

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @param i the row index
#' @param j the column index
#' @param ... other arguments (not used)
#' @export
setMethod(
    '[', c(x = "DistSum", i = "ANY", j = "ANY", drop = "missing"),
    function(x, i, j, ...) {
        .as.matrix(x)[i = i, j = j, ...]
    })

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @param i the array index
#' @param drop if TRUE, decrease the nb of dimensions when possible
#' @param ... other arguments (not used)
#' @export
setMethod(
    '[', c(x = "DistSum", i = "ANY", j = "missing", drop = "ANY"),
    function(x, i, j, ..., drop) {
        .as.matrix(x)[i = i, ..., drop = drop]
    })

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @param i the array index
#' @param ... other arguments (not used)
#' @export
setMethod(
    '[', c(x = "DistSum", i = "ANY", j = "missing", drop = "missing"),
    function(x, i, j, ...) {
        .as.matrix(x)[i = i, ...]
    })

.as.matrix <- function(x, whichFeatures = NULL) {
    stopifnot(inherits(x, "DistSum"))
    if (is.null(whichFeatures)) {
        return(x@pwFullDist)
    } else {
        if (length(unique(whichFeatures)) != length(whichFeatures)) {
            stop("whichFeatures array should not contain duplicates")
        }
        if (is.character(whichFeatures)) {
            whichFeatures <- vapply(
                X = whichFeatures,
                FUN.VALUE = integer(1),
                FUN = function(ft) {
                    ftInd <- which(featureNames(x) == ft)
                    if (length(ftInd) == 0){
                        stop("feature ", ft, " not found")
                    } 
                    ftInd
                }
            )
        } else if (is.numeric(whichFeatures)) {
            # check feature indices
            whichFeatures <- vapply(
                X = whichFeatures,
                FUN.VALUE = integer(1),
                FUN = function(ftInd) {
                    if (ftInd < 0 || ftInd > nFeatures(x)) {
                        stop("feature index out of bound [", ftInd,"]")
                    }
                    as.integer(ftInd)
                }
            )
        } else {
            stop("whichFeatures should be either NULL, a character vector, ",
                 "or a vector of indices")
        }
        if (length(whichFeatures) == nFeatures(x)) {
            return (x@pwFullDist)
        } else {
            dd <- Reduce(
                f = function(x, y){
                    x + y
                },
                x = x@pwDistPerDim[whichFeatures])
            
            return(dd)
        }
    }
}

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @param whichFeatures either an array of feature names, 
#' or an array of feature indices, or NULL
#' If NULL, the full distance (for all features) will be returned
#' If not NULL, `whichFeatures` array should not contain duplicates
#' @export
setMethod(
    "as.matrix", "DistSum",
    function(x, whichFeatures = NULL){
        return(.as.matrix(x, whichFeatures = whichFeatures))
    }
)

    




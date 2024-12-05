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
#' @slot pwDistPerFeature A `list` of `matrix` objects storing the contribution
#' of each feature (dimension) of the multidimensional distributions
#' to the full pairwise distance matrix. 
#' Note these matrices are not necessarily square symmetric matrices, 
#' as the `DistSum` could be occasionally used to store a given block of 
#' a bigger distance matrix. 
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
             pwDistPerFeature = "list"
         ),
         prototype = list(
             pwDistPerFeature = list()
         )
)

#' @importFrom methods is
.isValid <- function(object) {
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
    
    
    # check classes of slot list elements
    for (item in object@pwDistPerFeature) {
        if (!is(item, "matrix")) {
            return(paste0("pwDistPerFeature list contains other objects ",
                          "than `matrix` objects"))
        }
        if (!identical(dim(item), dim(object))) {
            return(paste0("pwDistPerFeature list contains matrices with ",
                          "different dimensions"))
        }
        if (!identical(dimnames(item), dimnames(object))) {
            return(paste0("pwDistPerFeature list contains matrices with ",
                          "different dim names"))
        }
    }
    return(TRUE)
}

setValidity("DistSum", function(object) {
    return(.isValid(object))
})

#' @rdname DistSum-class
#' @param object a `DistSum` object
#' @importMethodsFrom methods show
#' @importFrom methods validObject
setMethod(
    "show", "DistSum",
    function(object) {
        retVal <- methods::validObject(object)
        if (isTRUE(retVal)) {
            cat("`DistSum` object containing pairwise distances", 
                "between distributions\n", 
                "and their decomposition as a sum of feature contributions\n")
            cat("Matrix dimensions: ", dim(object), "\n")
            cat("Nb of features: ", nFeatures(object), "\n")
            cat("Feature names: ", featureNames(object), "\n")
            cat("Full distance matrix: \n")
            print(.as.matrix(object))
        } else {
            stop("Invalid object: ", retVal)
        }
    }
)

setGeneric("DistSum", function(object, ...) {
    standardGeneric("DistSum")
})

#' @rdname DistSum-class
#' @param object a `dist` object containing the full distance matrix
#' @importFrom methods validObject
#' @export
setMethod("DistSum",
          "matrix",
          function(object) {
              newObj <- methods::new("DistSum",
                                     pwDistPerFeature = list(object))


              if (isTRUE(methods::validObject(newObj))){
                  return(newObj)
              }
          })

#' @rdname DistSum-class
#' @param object a `list` of which each element is a `dist` object representing
#' the distance matrix for one feature.
#' @importFrom methods validObject
#' @export
setMethod(
    "DistSum", "list",
    function(object) {
        newObj <- methods::new("DistSum",
                               pwDistPerFeature = object)
        
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
        dd <- NULL
        if (length(x@pwDistPerFeature) > 0) {
            dd <- dim(x@pwDistPerFeature[[1]])    
        }
        return(dd)
    }
)

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @export
setMethod(
    "dimnames", "DistSum",
    function(x) {
        dn <- dimnames(x@pwDistPerFeature[[1]])
        return(dn)
    }
)

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @param value the new dimension names to be assigned
#' @importFrom methods validObject
#' @export
setMethod(
    "dimnames<-", c(x = "DistSum", value = "list"),
    function(x, value) {
        x@pwDistPerFeature <- lapply(
            x@pwDistPerFeature,
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
#' @importFrom methods validObject
#' @export
setMethod(
    "dimnames<-", c(x = "DistSum", value = "ANY"),
    function(x, value) {
        x@pwDistPerFeature <- lapply(
            x@pwDistPerFeature,
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
        nc <- ncol(x@pwDistPerFeature[[1]])
        return(nc)
    }
)

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @export
setMethod(
    "colnames", "DistSum",
    function(x) {
        cn <- colnames(x@pwDistPerFeature[[1]])
        return(cn)
    }
)

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @param value the new column names to be assigned
#' @importFrom methods validObject
#' @export
setMethod(
    "colnames<-", "DistSum",
    function(x, value) {
        x@pwDistPerFeature <- lapply(
            x@pwDistPerFeature,
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
        nr <- nrow(x@pwDistPerFeature[[1]])
        return(nr)
    }
)

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @export
setMethod(
    "rownames", "DistSum",
    function(x) {
        rn <- rownames(x@pwDistPerFeature[[1]])
        return(rn)
    }
)

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @param value the new column names to be assigned
#' @importFrom methods validObject
#' @export
setMethod(
    "rownames<-", "DistSum",
    function(x, value) {
        x@pwDistPerFeature <- lapply(
            x@pwDistPerFeature,
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
    return(length(x@pwDistPerFeature))
}

#' @rdname DistSum-class
#' @param object a `DistSum` object
#' @importMethodsFrom Biobase featureNames
#' @export
setMethod(
    "featureNames", "DistSum",
    function(object) {
     return(names(object@pwDistPerFeature))
    }
)

#' @rdname DistSum-class
#' @param object a `DistSum` object
#' @param value the new feature names to be assigned
#' @importMethodsFrom Biobase featureNames<-
#' @importFrom methods validObject
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
        names(object@pwDistPerFeature) <- value
        
        if (isTRUE(methods::validObject(object))) return(object)
    }
)

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @param i the row index
#' @param j the column index
#' @param ... other arguments (not used)
#' @param drop not supported (set to FALSE)
#' @importFrom methods validObject
#' @export
setMethod(
    '[', c(x = "DistSum", i = "ANY", j = "ANY", drop = "ANY"),
    function(x, i, j, ..., drop) {
        # note 'drop' should always set to false to keep slots as matrices and
        # list of matrices
        x@pwDistPerFeature <- lapply(
            x@pwDistPerFeature,
            FUN = function(mat){
                mat[i = i, j = j, ..., drop = FALSE]
            }
        )
        if (isTRUE(methods::validObject(x))) {
            return(x)    
        }
    })

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @param i the row index
#' @param j the column index
#' @param ... other arguments (not used)
#' @importFrom methods validObject
#' @export
setMethod(
    '[', c(x = "DistSum", i = "ANY", j = "ANY", drop = "missing"),
    function(x, i, j, ...) {
        # note 'drop' should always set to false to keep slots as matrices and
        # list of matrices
        x@pwDistPerFeature <- lapply(
            x@pwDistPerFeature,
            FUN = function(mat){
                mat[i = i, j = j, ..., drop = FALSE]
            }
        )
        if (isTRUE(methods::validObject(x))) {
            return(x)    
        }
    })

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @param i the array index
#' @param drop not supported (set to FALSE)
#' @param ... other arguments (not used)
#' @importFrom methods validObject
#' @export
setMethod(
    '[', c(x = "DistSum", i = "ANY", j = "missing", drop = "ANY"),
    function(x, i, j, ..., drop) {
        # note 'drop' should always set to false to keep slots as matrices and
        # list of matrices
        x@pwDistPerFeature <- lapply(
            x@pwDistPerFeature,
            FUN = function(mat){
                mat[i = i, ..., drop = FALSE]
            }
        )
        if (isTRUE(methods::validObject(x))) {
            return(x)    
        }
    })

#' @rdname DistSum-class
#' @param x a `DistSum` object
#' @param i the array index
#' @param ... other arguments (not used)
#' @importFrom methods validObject
#' @export
setMethod(
    '[', c(x = "DistSum", i = "ANY", j = "missing", drop = "missing"),
    function(x, i, j, ...) {
        # note 'drop' should always set to false to keep slots as matrices and
        # list of matrices
        x@pwDistPerFeature <- lapply(
            x@pwDistPerFeature,
            FUN = function(mat){
                mat[i = i, ..., drop = FALSE]
            }
        )
        if (isTRUE(methods::validObject(x))) {
            return(x)    
        }
        
    })

.as.matrix <- function(x, whichFeatures = NULL) {
    stopifnot(inherits(x, "DistSum"))
    if (is.null(whichFeatures)) {
        whichFeatures <- seq(nFeatures(x))
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
    }
    
    dd <- Reduce(
        f = function(x, y){
            x + y
        },
        x = x@pwDistPerFeature[whichFeatures])
    
    return(dd)    
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

    




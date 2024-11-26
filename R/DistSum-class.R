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
#' @rdname DistSum
#'
#' @description
#'
#' Class representing pairwise distances between multiple multidimensional 
#' distributions, when the distance is calculated as a sum of marginal 
#' distribution distances.
#' 
#' @slot pwFullDist An object of class `dist` storing the full pairwise
#' distance matrix. The `dist` object stores the triangular relevant part 
#' of the symmetrix, 0 diagonal, matrix only.
#' 
#' @slot pwDistPerDim A `list` of `dist` objects storing the contribution
#' of each dimension of the multidimensional distributions
#' to the full distance matrix
#'
#' @exportClass DistSum
#' 
#' @return nothing
#' @examples
#' 
#' nDistr <- 5
#' nFeat <- 7
#' 
#' M <- matrix(data = rnorm(nDistr * nFeat), ncol = nFeat)
#' rownames(M) <- paste0("distr", 1:nDistr)
#' colnames(M) <- paste0("feat", 1:nFeat)
#' DList <- lapply(colnames(M),
#'                 FUN = function(colName) {
#'                     dist(M[, colName, drop = FALSE])
#'                 })
#' 
#' names(DList) <- colnames(M)
#' D <- dist(M)
#' 
#' 
setClass("DistSum",
         slots = c(
             pwFullDist = "dist",
             pwDistPerDim = "list"
         ),
         prototype = list(
             pwFullDist = dist(numeric()),
             pwDistPerDim = list()
         )
)

# setGeneric("DistSum", function(object, ...) {
#     standardGeneric("DistSum")
# })
# 
# #' @rdname DistSum
# #' @param object a `dist` containing the full distance matrix
# #' @export
# #'
# setMethod(
#     "DistSum", "dist",
#     function(object) {
#         newObj <- methods::new("DistSum",
#                                pwFullDist = object,
#                                pwDistPerDim = list(object))
#         
#         
#         if (isTRUE(methods::validObject(newObj))){
#             return(newObj)
#         }             
#     }
# )
# 
# #' @rdname DistSum
# #' @param object a list of which each element is a `dist` object representing
# #' the distance matrix for one feature.
# #' @export
# #'
# setMethod(
#     "DistSum", "list",
#     function(object) {
#         pwFUllDist <- Reduce(
#             f = function(x, y){
#                 x + y
#             },
#             x = object)
#         
#         newObj <- methods::new("DistSum",
#                                pwFullDist = pwFUllDist,
#                                pwDistPerDim = object)
#                      
#         if (isTRUE(methods::validObject(newObj))){
#             return(newObj)
#         }
#     }
# )

setValidity("DistSum", function(object) {
    nDistr <- nDistr(object)
    dNames <- distrNames(object)
    if (!is.null(dNames)) {
        if (length(unique(dNames)) != 
            length(dNames)) {
            return("duplicate distribution names found")
        }
    }
    if (nDistr > 0) {
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
                                  "match distribution dimensions"))        
                }
            }
        } else {
            return("no features found in Dist object")
        }
    }
    
    #check classes of main objects
    if (!is(object@pwFullDist, "dist")) {
        return("pwFullDist object is not a `dist` object")
    }
    
    ret <- lapply(
        object@pwDistPerDim,
        FUN = function(item) {
            if (!is(item, "dist")) {
                return(paste0("pwDistPerDim list contains other objects ", 
                              "than `dist` objects"))
            }
       })
    
    return(TRUE)

})



#' @rdname DistSum
#' @param object a `DistSum` object
#'
#' @importMethodsFrom methods show
#'
#'
setMethod(
    "show", "DistSum",
    function(object) {
        cat("`DistSum` object containing pairwise distances", 
            "between distributions\n", 
            "and their decomposition as a sum of feature contributions\n",
            "Data:\n")
        cat("Nb of distributions: ", nDistr(object), "\n")
        if (nDistr(object) > 0) {
            cat("Nb of features: ", nFeatures(object), "\n")
            cat("Full distance matrix: \n")
            print(getPWDist(object))
        }
    }
)

#' @rdname DistSum
#' @param x a `DistSum` object
#' @export
nDistr <- function(x) {
    stopifnot(inherits(x, "DistSum"))
    return(attr(x@pwFullDist, "Size"))
}

#' @rdname DistSum
#' @param x a `DistSum` object
#' @export
nFeatures <- function(x) {
    stopifnot(inherits(x, "DistSum"))
    return(length(x@pwDistPerDim))
}

#' @rdname DistSum
#' @param x a `DistSum` object
#' @export
featureNames <- function(x) {
    stopifnot(inherits(x, "DistSum"))
    return(names(x@pwDistPerDim))
}

##' @rdname DistSum
##' @param x a `DistSum` object
##' @param value the new feature names to be assigned
##' @export
##'
"featureNames<-" <- function(x, value)
{
    stopifnot(inherits(x, "DistSum"))
    if (!is.character(value)) {
        stop("new feature names should be a character vector")
    }
    if (length(value) != nFeatures(x)) {
        stop("new feature names should have same length as nb of features")
    }
    names(x@pwDistPerDim) <- value
    
    if (isTRUE(methods::validObject(x))) return(x)
}

#' @rdname DistSum
#' @param x a `DistSum` object
#' @export
distrNames <- function(x) {
    stopifnot(inherits(x, "DistSum"))
    return(attr(x@pwFullDist, "Labels"))
}

##' @rdname DistSum
##' @param x a `DistSum` object
##' @param value the new feature names to be assigned
##' @export
##'
"distrNames<-" <- function(x, value)
{
    stopifnot(inherits(x, "DistSum"))
    if (!is.character(value)) {
        stop("new distribution names should be a character vector")
    }
    if (length(value) != nDistr(x)) {
        stop("new distribution names should have same length ",
             "as nb of distributions")
    }
    attr(x@pwFullDist, "Labels") <- value
    
    if (isTRUE(methods::validObject(x))) return(x)
}

#' @rdname DistSum
#' @param x a `DistSum` object
#' @param whichFeatures either an array of feature names, 
#' or an array of feature indices, or NULL
#' If NULL, the full disatnce (for all features) will be returned
#' @export
getPWDist <- function(x, whichFeatures = NULL) {
    stopifnot(inherits(x, "DistSum"))
    
    if (is.null(whichFeatures)) {
        return(x@pwFullDist)
    } else {
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
        if (length(whichFeatures) == nDistr(x)) {
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

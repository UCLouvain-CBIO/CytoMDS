# DistSum class

Class representing pairwise distances between multiple multidimensional
distributions, when the distance is calculated as a sum of marginal
distribution distances.

## Usage

``` r
# S4 method for class 'DistSum'
show(object)

# S4 method for class 'matrix'
DistSum(object)

# S4 method for class 'list'
DistSum(object)

# S4 method for class 'DistSum'
dim(x)

# S4 method for class 'DistSum'
dimnames(x)

# S4 method for class 'DistSum,list'
dimnames(x) <- value

# S4 method for class 'DistSum,ANY'
dimnames(x) <- value

# S4 method for class 'DistSum'
ncol(x)

# S4 method for class 'DistSum'
colnames(x)

# S4 method for class 'DistSum'
colnames(x) <- value

# S4 method for class 'DistSum'
nrow(x)

# S4 method for class 'DistSum'
rownames(x)

# S4 method for class 'DistSum'
rownames(x) <- value

nFeatures(x)

# S4 method for class 'DistSum'
featureNames(object)

# S4 method for class 'DistSum'
featureNames(object) <- value

# S4 method for class 'DistSum,ANY,ANY,ANY'
x[i, j, ..., drop = TRUE]

# S4 method for class 'DistSum,ANY,ANY,missing'
x[i, j, ..., drop = TRUE]

# S4 method for class 'DistSum,ANY,missing,ANY'
x[i, j, ..., drop = TRUE]

# S4 method for class 'DistSum,ANY,missing,missing'
x[i, j, ..., drop = TRUE]

# S4 method for class 'DistSum'
as.matrix(x, whichFeatures = NULL)

distByFeature(distObj)
```

## Arguments

- object:

  a `DistSum` object

- x:

  a `DistSum` object

- value:

  the new feature names to be assigned

- i:

  the array index

- j:

  the column index

- ...:

  other arguments (not used)

- drop:

  not supported (set to FALSE)

- whichFeatures:

  either an array of feature names, or an array of feature indices, or
  NULL If NULL, the full distance (for all features) will be returned If
  not NULL, `whichFeatures` array should not contain duplicates

- distObj:

  a `DistSum` object

## Value

nothing

a data.frame, with 3 columns:

- featureName : self explainatory

- distanceContrib : unidimensional distance along the corresponding
  feature

- percentage : percentage of feture distance w.r.t. full distance

## Slots

- `pwDistPerFeature`:

  A `list` of `matrix` objects storing the contribution of each feature
  (dimension) of the multidimensional distributions to the full pairwise
  distance matrix. Note these matrices are not necessarily square
  symmetric matrices, as the `DistSum` could be occasionally used to
  store a given block of a bigger distance matrix.

## Examples

``` r

# create a dummy distance matrix 
# to do this we use `nPoints` points 
# in an euclidian space of `nFeat` dimensions
nPoints <- 5
nFeat <- 7
M <- matrix(data = rnorm(nPoints * nFeat), ncol = nFeat)
rownames(M) <- paste0("point", 1:nPoints)
colnames(M) <- paste0("feat", 1:nFeat)

DList <- lapply(colnames(M),
FUN = function(colName) {
    D <- as.matrix(dist(
        M[, colName, drop = FALSE]))
    D
})

D <- Reduce(x = DList, f = function(A, B) A + B)

names(DList) <- colnames(M)

# Example of creating of a DistSum object based on the full distance matrix
distObj1 <- DistSum(D)
show(distObj1)
#> `DistSum` object containing pairwise distances between distributions
#>  and their decomposition as a sum of feature contributions
#> Matrix dimensions:  5 5 
#> Nb of features:  1 
#> Feature names:  
#> Full distance matrix: 
#>          point1    point2    point3    point4    point5
#> point1 0.000000  8.827115  7.094430  4.580522  8.285466
#> point2 8.827115  0.000000 10.959452 10.329413 10.876936
#> point3 7.094430 10.959452  0.000000  7.613967  6.125266
#> point4 4.580522 10.329413  7.613967  0.000000  9.235268
#> point5 8.285466 10.876936  6.125266  9.235268  0.000000

# Example of creation of a DistSum object based on a list of matrices
# representing the additive contribution of each feature
distObj2 <- DistSum(DList)

show(distObj2)
#> `DistSum` object containing pairwise distances between distributions
#>  and their decomposition as a sum of feature contributions
#> Matrix dimensions:  5 5 
#> Nb of features:  7 
#> Feature names:  feat1 feat2 feat3 feat4 feat5 feat6 feat7 
#> Full distance matrix: 
#>          point1    point2    point3    point4    point5
#> point1 0.000000  8.827115  7.094430  4.580522  8.285466
#> point2 8.827115  0.000000 10.959452 10.329413 10.876936
#> point3 7.094430 10.959452  0.000000  7.613967  6.125266
#> point4 4.580522 10.329413  7.613967  0.000000  9.235268
#> point5 8.285466 10.876936  6.125266  9.235268  0.000000

# getting dimensions
myDim <- dim(distObj2) # c(nPoints, nPoints)
ncols <- ncol(distObj2) # nPoints
nrows <- nrow(distObj2) # nPoints
nFeats <- nFeatures(distObj2) # nFeat
myFeatNames <- featureNames(distObj2) # paste0("feat", 1:nFeat)
myRowNames <- rownames(distObj2) # paste0("point", 1:nPoints)
myRowNames <- colnames(distObj2) # paste0("point", 1:nPoints)

# get full distance matrix
dd <- as.matrix(distObj2)

# get partial distance matrix for feature 1
dd1 <- as.matrix(distObj2, whichFeatures = 1)

# same thing, using feature name
dd1bis <- as.matrix(distObj2, whichFeatures = "feat1")

# getting partial distance for feature 1 & 2

ddPart <- as.matrix(distObj2, whichFeatures = colnames(M)[1:2])

# getting distance by feature
DF <- distByFeature(distObj2)
```

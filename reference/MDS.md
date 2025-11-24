# MDS class

Class representing Multi Dimensional Scaling (MDS) projection.

returns the value of the stress criterion, minimized by the SMACOF
algorithm.

returns a vector of nPoints dimension, containing the stress indicator
per point. The `stress` minimization criterion can indeed be allocated
per represented point. The more the stress of a particular point, the
less accurate its distances w.r.t. the other points.

## Usage

``` r
# S4 method for class 'MDS'
show(object)

nDim(x)

nPoints(x)

pwDist(x)

projections(x)

projDist(x)

stress(x)

spp(x)

eigenVals(x)

pctvar(x)

RSq(x)

RSqVec(x)

GoF(x)

smacofRes(x)
```

## Arguments

- object:

  a `MDS` object

- x:

  a `MDS` object

## Value

nothing

## Slots

- `nDim`:

  `numeric`, nb of dimensions of the projection

- `pwDist`:

  An object of class `dist` storing the triangular relevant part of the
  symmetric, zero diagonal pairwise distance matrix (nPoints \*
  nPoints), BEFORE projection.

- `proj`:

  The projection matrix, resulting from MDS

- `projDist`:

  An object of class `dist` storing the triangular relevant part of the
  symmetric, zero diagonal pairwise distance matrix (nPoints \*
  nPoints), AFTER projection.

- `eigen`:

  `numeric`, vector of `nDim` length, containing the eigen values of the
  PCA that is applied after the Smacof algorithm.

- `pctvar`:

  `numeric`, vector of `nDim` length, containing the percentage of
  explained variance per axis.

- `RSq`:

  `numeric`, vector of pseudo R square indicators, as a function of
  number of dimensions. `RSq[nDim]` is the global pseudo R square, as
  displayed on plots.

- `GoF`:

  `numeric`, vector of goodness of fit indicators, as a function of
  number of dimensions. `GoF[nDim]` is the global goodness of fit.

  Note pseudo R square and goodness of fit indicators are essentially
  the same indicator, only the definition of total sum of squares
  differ:

  - for pseudo RSq: TSS is calculated using the mean pairwise distance
    as minimum

  - for goodness of fit: TSS is calculated using 0 as minimum

- `smacofRes`:

  an object of class 'smacofB' containing the algorithmic optimization
  results, for example stress and stress per point, as returned by
  [`smacof::smacofSym()`](https://rdrr.io/pkg/smacof/man/smacofSym.html)
  method.

## Examples

``` r


nHD <- 10
nLD <- 2
nPoints <- 20 

# generate uniformly distributed points in 10 dimensions
points <- matrix(
    data = runif(n = nPoints * nHD),
    nrow = nPoints)
    
# calculate euclidian distances     
pwDist  <- dist(points)

# compute Metric MDS object by reaching a target pseudo RSquare
mdsObj <- computeMetricMDS(pwDist, targetPseudoRSq = 0.95)

show(mdsObj)
#> MDS object containing MDS projection (using Smacof algorithm)  data:
#> Nb of dimensions:  7 
#> Nb of points:  20 
#> Stress:  0.036138 
#> Pseudo RSquare:  0.958664 
#> Goodness of fit:  0.998694 

```

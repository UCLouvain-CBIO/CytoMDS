# metric MDS projection of sample

Multi-dimensional scaling projection of samples, using a distance matrix
as an input. The MDS algorithm is not the classical MDS (`cmdscale`
alike, aka Torgerson's algorithm), but is the SMACOF algorithm for
metric distances that are not necessarily euclidean. After having
obtained the projections on the `nDim` dimensions, we always apply svd
decomposition to visualize as first axes the ones that contain the most
variance of the projected dataset in `nDim` dimensions. Instead of being
provided directly by the user, the `nDim` parameter can otherwise be
found iteratively by finding the minimum `nDim` parameter that allows
the projection to reach a target pseudo RSquare. If this is the case,
the `maxDim` parameter is used to avoid looking for too big projection
spaces.

## Usage

``` r
computeMetricMDS(
  pwDist,
  whichChannels = NULL,
  nDim = NULL,
  seed = NULL,
  targetPseudoRSq = 0.95,
  maxDim = 128,
  ...
)
```

## Arguments

- pwDist:

  (`nSamples` rows, `nSamples` columns), previously calculated pairwise
  distances between samples, can be provided as :

  - a `DistSum` object

  - a `dist` object

  - a full symmetric square matrix, with 0. diagonal

- whichChannels:

  if `pwDist` has been provided as a `DistSum` object, a vector of
  channels to be included in the distances. In that case the distances
  have been computed as a sum of unidimensional distances for each
  channel, and the `DistSum` object allows to restrict the channel sets
  to be included in the distance accounting

- nDim:

  number of dimensions of projection, as input to SMACOF algorithm if
  not provided, will be found iteratively using `targetPseudoRSq`

- seed:

  seed to be set when launching SMACOF algorithm (e.g. when `init` is
  set to `"random"` but not only)

- targetPseudoRSq:

  target pseudo RSquare to be reached (only used when `nDim` is set to
  NULL)

- maxDim:

  in case `nDim` is found iteratively, maximum number of dimensions the
  search procedure is allowed to explore

- ...:

  additional parameters passed to SMACOF algorithm

## Value

an object of S4 class `MDS`

## Examples

``` r

library(CytoPipeline)

data(OMIP021Samples)

# estimate scale transformations 
# and transform the whole OMIP021Samples

transList <- estimateScaleTransforms(
    ff = OMIP021Samples[[1]],
    fluoMethod = "estimateLogicle",
    scatterMethod = "linearQuantile",
    scatterRefMarker = "BV785 - CD3")

OMIP021Trans <- CytoPipeline::applyScaleTransforms(
    OMIP021Samples, 
    transList)
    
# As there are only 2 samples in OMIP021Samples dataset,
# we create artificial samples that are random combinations of both samples

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

flowCore::pData(fsAll)$type <- factor(c("real", "real", rep("synthetic", 3)))
flowCore::pData(fsAll)$grpId <- factor(c("D1", "D2", rep("Agg", 3)))

# calculate all pairwise distances

pwDist <- pairwiseEMDDist(fsAll, 
                             channels = c("FSC-A", "SSC-A"),
                             verbose = FALSE)

# compute Metric MDS object with explicit number of dimensions
mdsObj <- computeMetricMDS(pwDist, nDim = 4, seed = 0)

dim <- nDim(mdsObj) # should be 4

#' # compute Metric MDS object by reaching a target pseudo RSquare
mdsObj2 <- computeMetricMDS(pwDist, seed = 0, targetPseudoRSq = 0.999)

```

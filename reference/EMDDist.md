# Calculate Earth Mover's distance between two samples

Calculate Earth Mover's distance between two samples

## Usage

``` r
EMDDist(
  x1,
  x2,
  channels = NULL,
  binSize = 0.05,
  minRange = -10,
  maxRange = 10,
  returnAll = FALSE
)
```

## Arguments

- x1:

  can be either a flowCore::flowFrame, or an expression matrix

- x2:

  can be either a flowCore::flowFrame, or an expression matrix

- channels:

  which channels (integer index(ices) or character(s)):

  - if it is a character vector, it can refer to either the channel
    names, or the marker names if x1 and x2 have been provided as
    flowCore::flowFrame

  - if it is a numeric vector, it refers to the indexes of channels in
    `x1`

  - if NULL : if `x1` and `x2` are provided as flowCore::flowFrames, all
    scatter and fluorescent channels of `x1` will be selected; if `x1`
    and `x2` are provided as expression matrices, all colnames of `x1`
    will be selected.

- binSize:

  size of equal bins to approximate the marginal distributions.

- minRange:

  minimum value taken when approximating the marginal distributions

- maxRange:

  maximum value taken when approximating the marginal distributions

- returnAll:

  If `TRUE`, distributions and marginal distribution distances are
  returned as well. Default = `FALSE`.

## Value

the Earth Mover's distance between `x1` and `x2`, which is calculated by
summing up all EMD approximates for the marginal distributions of each
channel

## Examples

``` r

library(CytoPipeline)
#> 
#> Attaching package: ‘CytoPipeline’
#> The following objects are masked from ‘package:Biobase’:
#> 
#>     pData, pData<-

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

# distance with itself (all channels at once)
# => should return 0
dist0 <- EMDDist(
    x1 = OMIP021Trans[[1]],
    x2 = OMIP021Trans[[1]])

# returning only distance, 2 channels
dist1 <- EMDDist(
    x1 = OMIP021Trans[[1]], 
    x2 = OMIP021Trans[[2]], 
    channels = c("FSC-A", "SSC-A"))

# using only one channel, passed by marker name
dist2 <- EMDDist(x1 = OMIP021Trans[[1]], 
                 x2 = OMIP021Trans[[2]], 
                 channels = c("BV785 - CD3"))

# using only one channel, passed by index
dist3 <- EMDDist(x1 = OMIP021Trans[[1]], 
                 x2 = OMIP021Trans[[2]], 
                 channels = 10)

dist2 == dist3
#> [1] TRUE
```

# Summary statistics per channel computation

Computation of summary statistic for selected channels, for all
flowFrames of a flowSet, or for all expression matrices of a list. This
method provides three different input modes:

- the user provides directly a flowCore::flowSet loaded in memory (RAM)

- the user provides directly a list of expression matrices of which the
  column names are the channel/marker names

- the user provides (1.) a number of samples `nSamples`; (2.) an ad-hoc
  function that takes as input an index between 1 and `nSamples`, and
  codes the method to load the corresponding expression matrix in
  memory;

## Usage

``` r
channelSummaryStats(
  x,
  loadExprMatrixFUN = NULL,
  loadExprMatrixFUNArgs = NULL,
  channels = NULL,
  statFUNs = stats::median,
  verbose = FALSE,
  BPPARAM = BiocParallel::SerialParam(),
  BPOPTIONS = BiocParallel::bpoptions(packages = c("flowCore"))
)
```

## Arguments

- x:

  can be:

  - a flowCore::flowSet

  - a list of expression matrices (Double matrix with named columns)

  - the number of samples (integer \>=1)

- loadExprMatrixFUN:

  the function used to translate an integer index into an expression
  matrix. In other words, the function should code how to load the
  `index`th expression matrix into memory. IMPORTANT: the expression
  matrix index should be the first function argument and should be named
  `exprMatrixIndex`.

- loadExprMatrixFUNArgs:

  (optional) a named list containing additional input parameters of
  `loadExprMatrixFUN()`

- channels:

  which channels needs to be included:

  - if it is a character vector, it can refer to either the channel
    names, or the marker names

  - if it is a numeric vector, it refers to the indices of channels in
    `fs`

  - if NULL, all scatter and fluorescent channels of `fs` \#' will be
    selected.

- statFUNs:

  a list (possibly of length one) of functions to call to calculate the
  statistics, or a simple function. This list can be named, in that
  case, these names will be transfered to the returned list.

- verbose:

  if `TRUE`, output a message after each single statistics calculation

- BPPARAM:

  sets the `BPPARAM` back-end to be used for the computation. If not
  provided, will use
  [`BiocParallel::SerialParam()`](https://rdrr.io/pkg/BiocParallel/man/SerialParam-class.html)
  (no task parallelization)

- BPOPTIONS:

  sets the BPOPTIONS to be passed to `bplapply()` function. Note that if
  you use a `SnowParams` back-end, you need to specify all the packages
  that need to be loaded for the different CytoProcessingStep to work
  properly (visibility of functions). As a minimum, the `flowCore`
  package needs to be loaded. (hence the default
  `BPOPTIONS = bpoptions(packages = c("flowCore"))` )

## Value

a list of named statistic matrices. In each stat matrix, the columns are
the channel statistics for all flowFrames of the flowSet. Exception: if
only one stat function (and not a list) is passed in `statFUNs`, the
return value is simplified to the stat matrix itself.

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

channelsOrMarkers <- c("FSC-A", "SSC-A", "BV785 - CD3")

# calculate mean for each 4 selected channels, for each 2 samples

channelMeans <- channelSummaryStats(
    OMIP021Trans,
    channels = channelsOrMarkers,
    statFUNs = mean)
    
# calculate median AND std deviation
# for each 4 selected channels, for each 2 samples

channelMedians <- channelSummaryStats(
    OMIP021Trans,
    channels = channelsOrMarkers,
    statFUNs = list("median" = stats::median, 
                    "std.dev" = stats::sd))
 
```

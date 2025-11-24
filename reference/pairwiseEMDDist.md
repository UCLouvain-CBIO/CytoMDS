# Pairwise Earth Mover's Distance calculation

Computation of all EMD between pairs of flowFrames belonging to a
flowSet. This method provides three different input modes:

- the user provides directly a flowCore::flowSet loaded in memory (RAM).

- the user provides directly a list of expression matrices loaded in
  RAM, of which the column names are the channel/marker names

- the user provides (1.) a number of samples `nSamples`; (2.) an ad-hoc
  function that takes as input an index between 1 and `nSamples`, and
  codes the method to load the corresponding expression matrix in
  memory; Optional row and column ranges can be provided to limit the
  calculation to a specific rectangle of the matrix. These i.e. can be
  specified as a way to split heavy calculations of large distance
  matrices on several computation nodes.

## Usage

``` r
pairwiseEMDDist(
  x,
  rowRange = c(1, nSamples),
  colRange = c(min(rowRange), nSamples),
  loadExprMatrixFUN = NULL,
  loadExprMatrixFUNArgs = NULL,
  channels = NULL,
  verbose = FALSE,
  BPPARAM = BiocParallel::SerialParam(),
  BPOPTIONS = BiocParallel::bpoptions(packages = c("flowCore")),
  binSize = 0.05,
  minRange = -10,
  maxRange = 10
)
```

## Arguments

- x:

  can be:

  - a flowCore::flowSet

  - a list of expression matrices (Double matrix with named columns)

  - the number of samples (integer \>=1)

- rowRange:

  the range of rows of the distance matrix to be calculated

- colRange:

  the range of columns of the distance matrix to be calculated

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

  which channels (integer index(ices) or character(s)):

  - if it is a character vector, it can refer to either the channel
    names, or the marker names

  - if it is a numeric vector, it refers to the indexes of channels in
    `fs`

  - if NULL all scatter and fluorescent channels of `fs` \#' will be
    selected

- verbose:

  if `TRUE`, output a message after each single distance calculation

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

- binSize:

  size of equal bins to approximate the marginal distributions.

- minRange:

  minimum value taken when approximating the marginal distributions

- maxRange:

  maximum value taken when approximating the marginal distributions

## Value

a distance matrix of pairwise distances (full symmetric with 0.
diagonal)

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

# calculate pairwise distances using only FSC-A & SSC-A channels
pwDist <- pairwiseEMDDist(
    x = OMIP021Trans,
    channels = c("FSC-A", "SSC-A"))
```

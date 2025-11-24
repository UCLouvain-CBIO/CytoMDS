# SampleMDS biplot wrapping

`ggplotSampleMDSWrapBiplots` calls `ggplotSampleMDS` repeatly to
generate biplots with different sets of external variables and align
them in a grid using the `patchwork` package, in a similar fashion as
[`ggplot2::facet_wrap()`](https://ggplot2.tidyverse.org/reference/facet_wrap.html)
does.

## Usage

``` r
ggplotSampleMDSWrapBiplots(
  mdsObj,
  extVariableList,
  ncol = NULL,
  nrow = NULL,
  byrow = NULL,
  displayLegend = TRUE,
  ...
)
```

## Arguments

- mdsObj:

  a MDS object, output of the
  [`computeMetricMDS()`](https://uclouvain-cbio.github.io/CytoMDS/reference/computeMetricMDS.md)
  method

- extVariableList:

  should be a named list of external variable matrices Each element of
  the list should be a matrix with named columns corresponding to the
  variables. The number of rows should be the same as the number of
  samples.

- ncol:

  passed to
  [`patchwork::wrap_plots()`](https://patchwork.data-imaginist.com/reference/wrap_plots.html)

- nrow:

  passed to
  [`patchwork::wrap_plots()`](https://patchwork.data-imaginist.com/reference/wrap_plots.html)

- byrow:

  passed to
  [`patchwork::wrap_plots()`](https://patchwork.data-imaginist.com/reference/wrap_plots.html)

- displayLegend:

  if FALSE, will de-active the legend display

- ...:

  additional parameters passed to
  [`ggplotSampleMDS()`](https://uclouvain-cbio.github.io/CytoMDS/reference/ggplotSampleMDS.md)
  (if used)

## Value

a ggplot object

## See also

[ggplotSampleMDS](https://uclouvain-cbio.github.io/CytoMDS/reference/ggplotSampleMDS.md),
[ggplotSampleMDSShepard](https://uclouvain-cbio.github.io/CytoMDS/reference/ggplotSampleMDSShepard.md),
[computeMetricMDS](https://uclouvain-cbio.github.io/CytoMDS/reference/computeMetricMDS.md)

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

# plot mds projection on axes 1 and 2,
# use 'group' for colour, 'type' for shape, and no label 

p_12 <- ggplotSampleMDS(
    mdsObj = mdsObj,
    pData = flowCore::pData(fsAll),
    projectionAxes = c(1,2),
    pDataForColour = "grpId",
    pDataForShape = "type")

# try to associate axes with median or std deviation of each channel
# => use bi-plots

extVarList <- channelSummaryStats(
    fsAll,
    channels = c("FSC-A", "SSC-A"),
    statFUNs = c("median" = stats::median, 
                 "std.dev" = stats::sd))

bpFull <- ggplotSampleMDSWrapBiplots(
    mdsObj = mdsObj,
    extVariableList = extVarList,
    pData = flowCore::pData(fsAll),
    projectionAxes = c(1,2),
    pDataForColour = "group",
    pDataForShape = "type",
    seed = 0)
```

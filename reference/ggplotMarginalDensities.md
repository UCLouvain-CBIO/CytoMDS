# Plot of channel intensity marginal densities

`ggplotMarginalDensities` uses ggplot2 to draw plots of marginal
densities of selected channels of a flowSet. If the flowSet contains
several flowFrames, events are concatenated together per group, or all
together in the absence of groups.

## Usage

``` r
ggplotMarginalDensities(
  x,
  sampleSubset,
  channels,
  pDataForColour,
  pDataForGroup,
  nEventInSubsample = Inf,
  seed = NULL,
  transList
)
```

## Arguments

- x:

  a
  [`flowCore::flowSet`](https://rdrr.io/pkg/flowCore/man/flowSet-class.html)
  (or a single
  [`flowCore::flowFrame`](https://rdrr.io/pkg/flowCore/man/flowFrame-class.html))

- sampleSubset:

  (optional) a logical vector, of the same length as `x`, indicating
  which flow frames to keep in the plot. Typically it is obtained
  through the evaluation of a logical condition the rows of
  `phenoData(fs)`.

- channels:

  (optional) - can be indices, or channel names, or markers.

- pDataForColour:

  (optional) which variable of `phenoData(fs)` will be used as colour
  aesthetic. Should be a character.

- pDataForGroup:

  (optional) which variable of `phenoData(fs)` will be used as group
  aesthetic. Should be a character. A separate marginal density will be
  calculated for each group and overlaid on the same channel density
  plots.

- nEventInSubsample:

  how many event to take (per flowFrame of the flowSet) for marginal
  density approximation.

- seed:

  if not null, used in subsampling.

- transList:

  a
  [`flowCore::transformList`](https://rdrr.io/pkg/flowCore/man/transformList-class.html)
  that will be applied before plotting.

## Value

a ggplot object

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

flowCore::pData(fsAll)$grpId <- factor(c("D1", "D2", rep("Agg", 3)))
flowCore::pData(fsAll)$lbl <- paste0("S", 1:5)

# plot densities, all samples together
p <- ggplotMarginalDensities(fsAll)

# plot densities, per sample
p <- ggplotMarginalDensities(fsAll, pDataForGroup = "lbl")

# plot densities, per sample and coloured by group
p <- ggplotMarginalDensities(
    fsAll, 
    pDataForGroup = "lbl",
    pDataForColour = "grpId")
```

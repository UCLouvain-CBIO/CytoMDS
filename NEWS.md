## CytoMDS 0.99

### CytoMDS 0.99.16
- added `lineWidth` parameter in `ggplotSampleMDSShepard()`
- running `plotly::ggplotly()` on `ggplotSampleMDSShepard()` output now 
displays row and column number for each distance point.

### CytoMDS 0.99.15
- corrected bug fix (error message) in `pwDist()` when `verbose=TRUE`

### CytoMDS 0.99.14
- re-factored code portions to replace, as much as possible, 
for loops by `apply()` family of functions.

### CytoMDS 0.99.13
- re-factored code portions to avoid growing lists incrementally

### CytoMDS 0.99.12
- removed `useBiocParallel` parameters from various stats functions 
(use BPPARAM = BiocParallel::SerialParam() as a default)
- implemented `MDS` class to store MDS projection results
- bi-plots now explicitly discard constant external variables (+warning) 
instead of raising an error without producing a plot
- implemented `ggplotMarginalDensities()`
- updated vignette with *Bodenmiller2012* dataset and more biological 
interpretation.

### CytoMDS 0.99.11
- re-factored package documentation file

### CytoMDS 0.99.10
- biplot now handles extVariables with missing values

### CytoMDS 0.99.9
- in `ggplotSampleMDS()` : add label layer after `geom_point()` (no more before)

### CytoMDS 0.99.8
- renamed `getChannelSummaryStats()` into `channelSummaryStats()`
- in `channelSummaryStats(), added support for `BiocParallel`, and allowed
for not loading the whole flowSet in memory at once.
- replaced NULL defaulted parameters with optional parameters
- added `displayPointLabels` argument to `ggplotSampleMDS()`
- added `displayLegend` argument to `ggplotSampleMDSWrapBiplots()`
- finalized creating vignette

### CytoMDS 0.99.7
- refactored the pairwise distance calculation code, by pre-computing the
unidimensional histograms and store them instead of recalculating them each
time a distance between 2 samples is calculated. This improves CPU time and
memory consumption.

### CytoMDS 0.99.6
- added `subset` argument in `ggplotSampleMDS()` and 
`ggplotSampleMDSWrapBiplots`

### CytoMDS 0.99.5
- renamed `getPairwiseEMDDist()` into `pairwiseEMDDist()`
- in `pairwiseEMDDist()`, added support for `BiocParallel`, and allowed
for not loading the whole flowSet in memory at once.

### CytoMDS 0.99.4
- in `getPairwiseEMDDist()`, added a second flowSet argument. When the two
flowSet arguments are non-null, distances are calculated for all sample pairs, 
where the first element comes from `fs`, 
and the second element comes from `fs2`.
- renamed `ggplotSamplesMDS` into `ggplotSampleMDS`
- renamed `ggplotSamplesMDSShepard` into `ggplotSampleMDSShepard`
- renamed `getChannelsSummaryStat` into `getChannelSummaryStats`
- new function `ggplotSampleMDSWrapBiplots()`

### CytoMDS 0.99.3
- new version of computeMetricMDS() which automatically sets 
the number of dimensions to reach a target *pseudo R squared*
- added ggplotly() functionality for output MDS plots
- in `ggplotSampleMDS()`, added `flipXAxis`, `flipYAxis` 
to possibly ease low dimensional projection comparisons
- in `ggplotSampleMDS()`, added `displayArrowLabels` to discard
the arrow labels in biplot. Also added `arrowThreshold`.
Moved arrow labels toward the end of the arrows.
- in `ggplotSampleMDS()` and `ggplotSampleMDSShepard()`: added 
`displayPseudoRSq` parameter.

### CytoMDS 0.99.2
- use global Rsquare as an indicator of quality of projection
- use %Var explained per axis

### CytoMDS 0.99.1
- in `ggplotSamplesMDS()`, added parameter `pDataForAdditionalLabelling`
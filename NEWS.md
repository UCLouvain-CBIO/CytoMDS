## CytoMDS 0.99.0

- Prior to Bioconductor submission

### CytoMDS 0.99.1
- in `ggplotSamplesMDS()`, added parameter `pDataForAdditionalLabelling`

### CytoMDS 0.99.2
- use global Rsquare as an indicator of quality of projection
- use %Var explained per axis

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

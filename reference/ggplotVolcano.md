# standard volcano plot

uses ggplot to draw a volcano plot

## Usage

``` r
ggplotVolcano(
  logFoldChanges,
  pValues,
  pointLabels = NULL,
  pointSizes = NULL,
  ggplotLabsArgs = NULL,
  pvalThresh = 0.05,
  pointColor = "blue",
  threshColor = "red",
  repelLabels = TRUE
)
```

## Arguments

- logFoldChanges:

  x axis values

- pValues:

  y axis values (will be -log10(pValues))

- pointLabels:

  labels to be attached to points

- pointSizes:

  point sizes

- ggplotLabsArgs:

  additional ggplot::labs() arguments (list)

- pvalThresh:

  p value threshold, if not NULL, a corresponding horizontal line will
  be plotted

- pointColor:

  color to use for points

- threshColor:

  color to use for the threshold line

- repelLabels:

  if TRUE uses
  [`ggrepel::geom_text_repel()`](https://ggrepel.slowkow.com/reference/geom_text_repel.html)

## Value

a ggplot object

## Examples

``` r

LFC <- c(-3,-2.7,-2.6,-1.5,-0.3,0.1,0.7,1.4,1.9,2.6)
pval <- c(0.004, 0.03, 0.022, 0.06, 0.4, 0.7, 0.3, 0.055, 0.045, 0.02)
labels <- paste0("p", (1:10)[sample(1:10, 10)])

set.seed(1)

# ggplotVolcano with default params
p <- ggplotVolcano(LFC, pval, labels)

# ggplotVolcano with no thresh
p <- ggplotVolcano(LFC, pval, labels, pvalThresh = NULL)

# ggplotVolcano with other thresh
p <- ggplotVolcano(LFC, pval, labels, pvalThresh = 0.01)

# ggplotVolcano with other ggplot_labs
p <- ggplotVolcano(LFC, pval, labels,
    ggplotLabsArgs = list(x = "my log2 fold",
                          y = "my log10 pval",
                          title = "My wonderful volcano!"))
                          
# ggplotVolcano with colors
p <- ggplotVolcano(LFC, pval, labels,
                   pointColor = "green", threshColor = "purple")
                   
# ggplotVolcano with point sizes
logNCells <- c(0.5, 1, 0.8, 2, 1.6, 1.2, 2, 1.2, 0.4, 0.2)
ggplotLabsArgs <- list(size = "nb cells (log10)")
p <- ggplotVolcano(LFC, pval, labels,
                   pointSizes = logNCells,
                   ggplotLabsArgs = ggplotLabsArgs)
```

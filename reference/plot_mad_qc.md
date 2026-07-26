# Plot an existing MAD QC report

Plot an existing MAD QC report

## Usage

``` r
plot_mad_qc(
  qc,
  metrics = NULL,
  type = c("distribution", "index"),
  facet_by = NULL,
  show_thresholds = TRUE,
  show_legend = TRUE
)
```

## Arguments

- qc:

  A tidy report created by
  [`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md).
  Calculation groups are read from report attributes when group-specific
  thresholds must be displayed.

- metrics:

  `NULL` for all metrics or a character vector to plot.

- type:

  A distribution or observation-index view.

- facet_by:

  `NULL` or metadata columns used only to arrange the plot. Faceting
  never changes MAD calculations or thresholds.

- show_thresholds:

  Show active lower and upper thresholds?

- show_legend:

  Show legends?

## Value

A ggplot object. The plot is returned without being printed.

## Details

Missing QC values are shown as blue triangles at an artificial panel
floor. The report itself remains unchanged. Lower thresholds are dashed
and upper thresholds are dot-dashed. No thresholds are recalculated.

## Examples

``` r
d <- data.frame(library_size = c(rnorm(19, 1e6, 1e5), NA))
qc <- mad_qc(d, c(library_size = "lower"))
plot_mad_qc(qc)
```

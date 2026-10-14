# Plot the relationship between two QC metrics

Show paired raw values and per-batch Spearman correlations. The
correlation and finite-pair count are descriptive; no p-values are
calculated. Requires `ggplot2`.

## Usage

``` r
plot_mad_qc_scatter(data, x, y, batch = NULL)
```

## Arguments

- data:

  Annotated data frame or Seurat object returned by
  [`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md).

- x, y:

  Character names of two distinct annotated numeric metrics.

- batch:

  Optional non-missing vector with one batch label per observation. If
  omitted, plot all observations in one panel.

## Value

A `ggplot` object with one panel per batch.

## Examples

``` r
if (requireNamespace("ggplot2", quietly = TRUE)) {
  x <- data.frame(counts = c(1:5, 20), pct_mito = c(2:6, 30))
  x <- mad_qc(x, c(counts = "lower", pct_mito = "upper"), verbose = FALSE)
  plot_mad_qc_scatter(x, "counts", "pct_mito")
}
```

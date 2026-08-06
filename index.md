# veryMAD

![veryMAD R package logo](reference/figures/veryMAD-logo.svg)

veryMAD is a small R package for explicit MAD-based quality-control
annotation and robust MAD scaling. It calculates flags and thresholds
but never filters observations or infers biological causes.

## Installation

``` r

pak::pak("Thokas99/veryMAD")
```

## MAD QC

Rows are observations. Select QC columns and their directions
explicitly.

``` r

library(veryMAD)

metadata <- data.frame(
  library_size = c(2e6, 31e6, 32e6, 30e6, 33e6),
  mapping_rate = c(0.50, 0.92, 0.94, 0.93, 0.95),
  row.names = paste0("sample", 1:5)
)

annotated <- mad_qc(
  metadata,
  metrics = c(library_size = "lower", mapping_rate = "lower"),
  transform = c(library_size = "log1p")
)
annotated
```

For thresholds and one-row-per-observation flags:

``` r

report <- mad_qc(
  metadata,
  metrics = c(library_size = "lower", mapping_rate = "lower"),
  transform = c(library_size = "log1p"),
  output = "report",
  verbose = FALSE
)
report$thresholds
report$flags
```

Use [`split()`](https://rdrr.io/r/base/split.html) and separate calls
when stratified QC is scientifically needed. `min_n` is a computational
safeguard, not a biological rule.

## MAD scaling

``` r

mad_scale(matrix_data, margin = 1) # row-wise
mad_z_score(matrix_data, margin = 2) # column-wise alias
```

`margin = 1` scales rows; `margin = 2` scales columns. Vectors are
scaled directly, and data frames are returned as matrices.

## Interpretation

MAD thresholds are adaptive statistical heuristics. Flags are not
diagnoses, upper-tail count flags are not doublet calls, and veryMAD
does not filter data automatically. Users choose metrics, directions,
and transformations.

## Documentation

See the [reference](https://thokas99.github.io/veryMAD/reference/) and
[guides](https://thokas99.github.io/veryMAD/articles/).

# Explicit MAD quality-control annotation

Rows are observations. Select numeric metric columns and a tail
direction for each. veryMAD estimates one threshold per metric across
the input rows and returns annotated metadata or a `verymad_qc` report.

## Usage

``` r
mad_qc(
  data,
  metrics,
  nmads = 3,
  transform = "none",
  output = c("annotate", "report"),
  min_n = 5,
  zero_mad = c("na", "zero", "error"),
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- data:

  A data frame, tibble, numeric matrix, or Seurat object.

- metrics:

  A non-empty named character vector with unique metric names; values
  are `"lower"`, `"upper"`, or `"both"`.

- nmads:

  Positive number of MADs from the median for the limits.

- transform:

  A scalar transformation or named partial overrides. Supported values
  are `"none"`, `"log1p"`, and `"log10"`.

- output:

  Return annotated data by default; use `"report"` for flags,
  thresholds, and settings.

- min_n:

  Minimum number of finite, non-missing observations required.

- zero_mad:

  Handling of a zero MAD: `"na"`, `"zero"`, or `"error"`.

- overwrite:

  Replace existing flag columns?

- verbose:

  Print a summary?

## Value

Annotated input or a `verymad_qc` list with `flags`, `thresholds`, and
`settings`. Threshold reports include finite observations used and
per-metric flagged and unresolved counts.

## Examples

``` r
qc_metadata <- data.frame(
  sample = paste0("sample_", 1:6),
  library_size = c(2.4e6, 2.5e6, 2.6e6, 2.7e6, 0.5e6, 2.5e6),
  pct_mito = c(4, 5, 6, 5, 7, 30)
)
annotated <- mad_qc(
  qc_metadata,
  metrics = c(library_size = "lower", pct_mito = "upper"),
  verbose = FALSE
)
annotated[, c("sample", "mad_qc_outlier")]
#>     sample mad_qc_outlier
#> 1 sample_1          FALSE
#> 2 sample_2          FALSE
#> 3 sample_3          FALSE
#> 4 sample_4          FALSE
#> 5 sample_5           TRUE
#> 6 sample_6           TRUE

report <- mad_qc(
  qc_metadata,
  metrics = c(library_size = "lower", pct_mito = "upper"),
  output = "report",
  verbose = FALSE
)
report$thresholds[, c("metric", "n_used", "n_flagged", "n_unresolved")]
#>         metric n_used n_flagged n_unresolved
#> 1 library_size      6         1            0
#> 2     pct_mito      6         1            0
```

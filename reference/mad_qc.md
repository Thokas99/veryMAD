# Annotate MAD-based quality-control flags

Rows are observations. For each selected metric, estimate limits around
the median using `nmads` scaled MADs, then flag the selected tail.
Transform metrics before estimation, fit limits by `batch`, or select
reference rows with `subset`. Return flags for every observation.

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
  verbose = TRUE,
  batch = NULL,
  subset = NULL
)
```

## Arguments

- data:

  A data frame, numeric matrix, or Seurat object. For Seurat, metrics
  must already be present in cell metadata.

- metrics:

  A non-empty named character vector with unique numeric column names
  and values `"lower"`, `"upper"`, or `"both"`.

- nmads:

  Positive finite number of scaled MADs from the median used for each
  limit. The MAD uses the normal-consistency constant 1.4826.

- transform:

  One method for all metrics or named overrides; unspecified metrics use
  `"none"`. Supported methods are `"none"`, `"log1p"`, and `"log10"`.

- output:

  `"annotate"` adds flags to the input; `"report"` returns flags,
  limits, and settings without modifying the input.

- min_n:

  Positive integer minimum number of non-missing reference observations
  per metric and batch group (default 5). Smaller groups have unresolved
  flags.

- zero_mad:

  Policy when a group's MAD is zero: `"na"` (default) leaves flags
  unresolved, `"zero"` sets the limit at the median, and `"error"`
  stops.

- overwrite:

  Replace existing QC flag columns when annotating?

- verbose:

  Print a summary?

- batch:

  Optional non-missing vector with one label per observation. Estimate
  separate limits within each label.

- subset:

  Optional non-missing logical vector selecting reference observations
  for estimation. Apply the resulting limits to every row.

## Value

With `output = "annotate"`, returns the input with per-metric and
combined flags; matrix input returns a data frame and Seurat input
remains a Seurat object. With `output = "report"`, returns a
`verymad_qc` list containing `flags`, `thresholds`, and `settings`.
Threshold rows report limits on transformed and raw scales, status, and
used, flagged, and unresolved counts. A batch adds one row per metric
and batch.

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

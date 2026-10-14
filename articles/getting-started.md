# Getting started with veryMAD

[`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md)
estimates one threshold per metric from the median and scaled MAD.
Choose the metrics, tail directions, and transformations. Annotation
retains all rows; report mode returns flags, thresholds, and settings
without modifying the input.

> **Note**
>
> A flag identifies an outlier in the supplied distribution. Review each
> metric in context before deciding whether to filter an observation.

## Annotate metrics

The example uses percentage points for `pct_mito`: 8 denotes 8%.

``` r

sample_metadata <- data.frame(
  sample = paste0("sample_", 1:8),
  library_size = c(2.4e6, 2.5e6, 2.6e6, 2.7e6, 2.8e6, 2.6e6, 0.45e6, 2.5e6),
  detected_genes = c(12000, 12500, 13100, 12800, 13400, 12600, 1500, 12300),
  pct_mito = c(4, 5, 6, 5, 7, 6, 8, 29)
)

qc_directions <- c(
  library_size = "lower",
  detected_genes = "lower",
  pct_mito = "upper"
)

qc_annotated <- mad_qc(
  sample_metadata,
  metrics = qc_directions,
  transform = c(library_size = "log1p", detected_genes = "log1p"),
  verbose = FALSE
)

qc_annotated[, c("sample", "library_size_mad_outlier",
                 "detected_genes_mad_outlier", "pct_mito_mad_outlier",
                 "mad_qc_outlier")]
#>     sample library_size_mad_outlier detected_genes_mad_outlier
#> 1 sample_1                    FALSE                      FALSE
#> 2 sample_2                    FALSE                      FALSE
#> 3 sample_3                    FALSE                      FALSE
#> 4 sample_4                    FALSE                      FALSE
#> 5 sample_5                    FALSE                      FALSE
#> 6 sample_6                    FALSE                      FALSE
#> 7 sample_7                     TRUE                       TRUE
#> 8 sample_8                    FALSE                      FALSE
#>   pct_mito_mad_outlier mad_qc_outlier
#> 1                FALSE          FALSE
#> 2                FALSE          FALSE
#> 3                FALSE          FALSE
#> 4                FALSE          FALSE
#> 5                FALSE          FALSE
#> 6                FALSE          FALSE
#> 7                FALSE           TRUE
#> 8                 TRUE           TRUE
```

Each metric flag is `TRUE`, `FALSE`, or `NA` when its threshold is
unresolved. The combined flag is `TRUE` if any metric is flagged;
otherwise it is `NA` if any metric remains unresolved.

## Inspect thresholds

Use `output = "report"` to inspect thresholds and flags without adding
columns to the input.

``` r

qc_report <- mad_qc(
  sample_metadata,
  metrics = qc_directions,
  transform = c(library_size = "log1p", detected_genes = "log1p"),
  output = "report",
  verbose = FALSE
)

qc_report$thresholds[, c("metric", "direction", "lower_raw",
                         "upper_raw", "status", "n_used", "n_flagged",
                         "n_unresolved")]
#>           metric direction  lower_raw upper_raw status n_used n_flagged
#> 1   library_size     lower 2148451.65        NA     ok      8         1
#> 2 detected_genes     lower   10908.88        NA     ok      8         1
#> 3       pct_mito     upper         NA   10.4478     ok      8         1
#>   n_unresolved
#> 1            0
#> 2            0
#> 3            0
```

Threshold rows report limits on the transformed and raw scales. `n_used`
is the number of non-missing reference observations; `min_n` defaults to
5.

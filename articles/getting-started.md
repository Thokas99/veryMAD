# Getting started with veryMAD

`veryMAD` annotates observations using robust, median absolute deviation
(MAD) thresholds. You choose the metrics, the direction of a poor value,
and any transformation. The input rows stay in the result.

## A basic QC workflow

The example uses a small sequencing metadata table. Low library size and
low gene detection suggest poor complexity; high mitochondrial content
is also a common QC concern.

``` r

sample_metadata <- data.frame(
  sample = paste0("sample_", 1:8),
  library_size = c(2.4e6, 2.5e6, 2.6e6, 2.7e6, 2.8e6, 2.6e6, 0.45e6, 2.5e6),
  detected_genes = c(12000, 12500, 13100, 12800, 13400, 12600, 1500, 12300),
  pct_mito = c(.04, .05, .06, .05, .07, .06, .08, .29)
)

# Describe which tail is a potential QC problem for each metric
qc_directions <- c(
  library_size = "lower",
  detected_genes = "lower",
  pct_mito = "upper"
)

# Add one flag per metric and a combined flag to the metadata
qc_annotated <- mad_qc(
  sample_metadata,
  metrics = qc_directions,
  verbose = FALSE
)
```

The result keeps all eight samples and adds logical flag columns.
Inspect the columns that matter for a quick review:

``` r

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

`mad_qc_outlier` is `TRUE` when at least one selected metric is flagged.
veryMAD does not remove those rows, so the decision about filtering or
follow-up stays with the analysis workflow.

## Inspect thresholds with a report

Use `output = "report"` when you want thresholds and flags in a compact
object instead of adding columns to the input table.

``` r

qc_report <- mad_qc(
  sample_metadata,
  metrics = qc_directions,
  output = "report",
  verbose = FALSE
)

qc_report$thresholds[, c("metric", "direction", "lower_raw",
                         "upper_raw", "status")]
#>           metric direction  lower_raw upper_raw status
#> 1   library_size     lower 2105220.00        NA     ok
#> 2 detected_genes     lower   10770.88        NA     ok
#> 3       pct_mito     upper         NA  0.104478     ok
qc_report$flags[, c("id", "mad_qc_outlier")]
#>   id mad_qc_outlier
#> 1  1          FALSE
#> 2  2          FALSE
#> 3  3          FALSE
#> 4  4          FALSE
#> 5  5          FALSE
#> 6  6          FALSE
#> 7  7           TRUE
#> 8  8           TRUE
```

The threshold table records the calculation status for each metric.
Missing or undersized metrics remain visible there instead of being
silently turned into a filtering decision.

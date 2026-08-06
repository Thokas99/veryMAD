# Getting started with veryMAD

`veryMAD` accepts observation metadata and requires users to select
numeric metrics, directions, and transformations explicitly. This
example uses a small synthetic bulk RNA-seq cohort with 36 libraries.
The values are only illustrative, but the workflow is the same for a
real metadata table.

``` r

set.seed(101)
n_libraries <- 36L
bulk_metadata <- data.frame(
  sample_id = sprintf("bulk_%02d", seq_len(n_libraries)),
  library_size = round(rlnorm(n_libraries, log(2.5e6), 0.18)),
  mapping_rate = pmin(.995, pmax(.70, rnorm(n_libraries, .94, .015))),
  duplication_rate = pmin(.95, pmax(.10, rnorm(n_libraries, .35, .04)))
)
bulk_metadata$library_size[c(2, 11)] <- c(180000, 220000)
bulk_metadata$mapping_rate[c(5, 24)] <- c(.70, .73)
bulk_metadata$duplication_rate[31] <- .82
rownames(bulk_metadata) <- bulk_metadata$sample_id

bulk_annotated <- mad_qc(
  bulk_metadata,
  metrics = c(
    library_size = "lower",
    mapping_rate = "lower",
    duplication_rate = "upper"
  ),
  transform = c(library_size = "log1p"),
  min_n = 10,
  verbose = FALSE
)

bulk_annotated[bulk_annotated$mad_qc_outlier %in% TRUE,
  c("sample_id", "library_size_mad_outlier", "mapping_rate_mad_outlier",
    "duplication_rate_mad_outlier", "mad_qc_outlier")]
#>         sample_id library_size_mad_outlier mapping_rate_mad_outlier
#> bulk_02   bulk_02                     TRUE                    FALSE
#> bulk_05   bulk_05                    FALSE                     TRUE
#> bulk_11   bulk_11                     TRUE                    FALSE
#> bulk_24   bulk_24                    FALSE                     TRUE
#> bulk_31   bulk_31                    FALSE                    FALSE
#>         duplication_rate_mad_outlier mad_qc_outlier
#> bulk_02                        FALSE           TRUE
#> bulk_05                        FALSE           TRUE
#> bulk_11                        FALSE           TRUE
#> bulk_24                        FALSE           TRUE
#> bulk_31                         TRUE           TRUE
```

The annotation preserves all 36 observations and adds one logical column
per metric, plus `mad_qc_outlier`. The combined flag is `TRUE` when any
selected metric is flagged; the input rows are never filtered.

Use `output = "report"` when thresholds and observation IDs should be
kept in a compact object rather than added to the metadata table:

``` r

bulk_report <- mad_qc(
  bulk_metadata,
  metrics = c(
    library_size = "lower",
    mapping_rate = "lower",
    duplication_rate = "upper"
  ),
  transform = c(library_size = "log1p"),
  output = "report",
  min_n = 10,
  verbose = FALSE
)

bulk_report$thresholds[, c("metric", "direction", "transform", "lower_raw",
                           "upper_raw", "status")]
#>             metric direction transform    lower_raw upper_raw status
#> 1     library_size     lower     log1p 1.431377e+06        NA     ok
#> 2     mapping_rate     lower      none 8.961747e-01        NA     ok
#> 3 duplication_rate     upper      none           NA 0.4744178     ok
head(bulk_report$flags)
#>        id library_size mapping_rate duplication_rate mad_qc_outlier
#> 1 bulk_01        FALSE        FALSE            FALSE          FALSE
#> 2 bulk_02         TRUE        FALSE            FALSE           TRUE
#> 3 bulk_03        FALSE        FALSE            FALSE          FALSE
#> 4 bulk_04        FALSE        FALSE            FALSE          FALSE
#> 5 bulk_05        FALSE         TRUE            FALSE           TRUE
#> 6 bulk_06        FALSE        FALSE            FALSE          FALSE
```

# Interpreting MAD QC flags

Choose `lower` for unusually small values, `upper` for unusually large
values, and `both` for either tail. Count-like metrics may benefit from
explicit `log1p`; bounded rates are often most interpretable on their
raw scale.

The following synthetic single-cell metadata contains 180 cells and four
common QC metrics. A few cells are deliberately made unusual so that the
resulting annotation is easy to inspect.

``` r

set.seed(202)
n_cells <- 180L
sc_metadata <- data.frame(
  cell_id = sprintf("cell_%03d", seq_len(n_cells)),
  nCount_RNA = round(rlnorm(n_cells, log(4000), .45)),
  nFeature_RNA = round(rlnorm(n_cells, log(1200), .30)),
  percent_mito = pmin(.40, pmax(.005, rbeta(n_cells, 2, 18))),
  doublet_score = pmin(.99, pmax(.01, rnorm(n_cells, .12, .04)))
)
sc_metadata$nCount_RNA[c(1, 2)] <- c(140, 180)
sc_metadata$nFeature_RNA[c(1, 2)] <- c(85, 120)
sc_metadata$percent_mito[c(3, 4)] <- c(.62, .55)
sc_metadata$doublet_score[5] <- .98
rownames(sc_metadata) <- sc_metadata$cell_id

sc_annotated <- mad_qc(
  sc_metadata,
  metrics = c(
    nCount_RNA = "lower",
    nFeature_RNA = "lower",
    percent_mito = "upper",
    doublet_score = "upper"
  ),
  transform = c(nCount_RNA = "log1p", nFeature_RNA = "log1p"),
  min_n = 20,
  verbose = FALSE
)

sc_annotated[sc_annotated$mad_qc_outlier %in% TRUE,
  c("cell_id", "nCount_RNA_mad_outlier", "nFeature_RNA_mad_outlier",
    "percent_mito_mad_outlier", "doublet_score_mad_outlier",
    "mad_qc_outlier")]
#>           cell_id nCount_RNA_mad_outlier nFeature_RNA_mad_outlier
#> cell_001 cell_001                   TRUE                     TRUE
#> cell_002 cell_002                   TRUE                     TRUE
#> cell_003 cell_003                  FALSE                    FALSE
#> cell_004 cell_004                  FALSE                    FALSE
#> cell_005 cell_005                  FALSE                    FALSE
#> cell_050 cell_050                  FALSE                    FALSE
#> cell_073 cell_073                  FALSE                    FALSE
#>          percent_mito_mad_outlier doublet_score_mad_outlier mad_qc_outlier
#> cell_001                    FALSE                     FALSE           TRUE
#> cell_002                    FALSE                     FALSE           TRUE
#> cell_003                     TRUE                     FALSE           TRUE
#> cell_004                     TRUE                     FALSE           TRUE
#> cell_005                    FALSE                      TRUE           TRUE
#> cell_050                     TRUE                     FALSE           TRUE
#> cell_073                     TRUE                     FALSE           TRUE
```

This is an annotation step, not a filtering step. The full table remains
available in `sc_annotated`; a downstream workflow can decide how to
review or use the flagged cells.

``` r

sc_report <- mad_qc(
  sc_metadata,
  metrics = c(
    nCount_RNA = "lower",
    nFeature_RNA = "lower",
    percent_mito = "upper",
    doublet_score = "upper"
  ),
  transform = c(nCount_RNA = "log1p", nFeature_RNA = "log1p"),
  output = "report",
  min_n = 20,
  verbose = FALSE
)

sc_report$thresholds[, c("metric", "direction", "transform", "median", "mad",
                         "lower_raw", "upper_raw", "status")]
#>          metric direction transform     median        mad lower_raw upper_raw
#> 1    nCount_RNA     lower     log1p 8.30511100 0.47135824  982.4200        NA
#> 2  nFeature_RNA     lower     log1p 7.07199691 0.31699966  454.3188        NA
#> 3  percent_mito     upper      none 0.07467968 0.05235087        NA 0.2317323
#> 4 doublet_score     upper      none 0.12177996 0.04920403        NA 0.2693921
#>   status
#> 1     ok
#> 2     ok
#> 3     ok
#> 4     ok
```

`min_n` is a computational safeguard. A zero MAD is reported explicitly
and can be handled as `"na"`, `"zero"`, or `"error"`. These flags are
statistical heuristics, not biological diagnoses or automatic filtering
decisions. The overall `mad_qc_outlier` flag is always returned as a
compact any-selected- metric summary; it is `NA` when no metric is
flagged but at least one metric could not be evaluated.

For stratified analyses, split the metadata explicitly and call
[`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md) on
each subset. veryMAD does not infer batches, conditions, clusters, or
causes.

# Explicit observation-level QC with veryMAD

## Data model

Rows are observations and columns are metadata or candidate QC metrics.
Every metric and direction is user-selected; veryMAD never selects
numeric columns, recognizes naming conventions, or filters observations.
The compact report is a list containing `flags`, `thresholds`, and
`settings`.

## Bulk observations

``` r

bulk <- data.frame(
  library_size = c(1e6, 30e6, 32e6, 31e6, 33e6),
  mapping_rate = c(.50, .92, .94, .93, .95),
  batch = c("a", "a", "a", "b", "b"),
  row.names = paste0("sample", 1:5)
)
metrics <- c(library_size = "lower", mapping_rate = "lower")
flagged <- mad_qc_bulk(bulk, metrics,
  transform = c(mapping_rate = "none"), verbose = FALSE)
report <- attr(flagged, "mad_qc")
report$thresholds
#> # A tibble: 2 × 14
#>   metric  direction transform n_observations n_evaluated n_missing median    mad
#>   <chr>   <chr>     <chr>              <int>       <int>     <int>  <dbl>  <dbl>
#> 1 librar… lower     log1p                  5           5         0  17.2  0.0486
#> 2 mappin… lower     none                   5           5         0   0.93 0.0148
#> # ℹ 6 more variables: lower <dbl>, upper <dbl>, lower_raw <dbl>,
#> #   upper_raw <dbl>, n_outliers <int>, outlier_proportion <dbl>
```

The same call accepts a matrix whose rows are observations. A
gene-by-sample expression matrix is not automatically a QC table.

``` r

mad_qc_bulk(as.matrix(bulk[c("library_size", "mapping_rate")]), metrics,
  transform = c(mapping_rate = "none"), verbose = FALSE)
#>         library_size mapping_rate library_size_mad_outlier
#> sample1      1.0e+06         0.50                     TRUE
#> sample2      3.0e+07         0.92                    FALSE
#> sample3      3.2e+07         0.94                    FALSE
#> sample4      3.1e+07         0.93                    FALSE
#> sample5      3.3e+07         0.95                    FALSE
#>         mapping_rate_mad_outlier
#> sample1                     TRUE
#> sample2                    FALSE
#> sample3                    FALSE
#> sample4                    FALSE
#> sample5                    FALSE
```

## Single-cell metadata

``` r

cells <- data.frame(
  nCount_RNA = c(NA, 500, 520, 540, 5000),
  nFeature_RNA = c(NA, 200, 210, 220, 1000),
  percent.mt = c(NA, 3, 4, 5, 30),
  row.names = paste0("cell", 1:5)
)
cell_metrics <- c(nCount_RNA = "lower", nFeature_RNA = "lower", percent.mt = "upper")
sc_report <- mad_qc_sc(cells, cell_metrics,
  transform = c(percent.mt = "none"), action = "report", verbose = FALSE)
annotated <- mad_qc_sc(cells, cell_metrics,
  transform = c(percent.mt = "none"), verbose = FALSE)
annotated$mad_qc_outlier
#> [1]    NA FALSE FALSE FALSE  TRUE
```

`TRUE` means any evaluated metric failed; `FALSE` means all were
evaluated and passed; `NA` means none failed but at least one was
unavailable. Individual flags stay in `sc_report$flags`. With a Seurat
object, the same annotation call adds only `mad_qc_outlier` through
optional `SeuratObject` and changes no cells, assays, reductions, or
identities.

## Transformations and groups

The wrapper default is `log1p` for every selected metric. Named values
are partial overrides, so bounded rates can remain on their original
scale without repeating every count metric. Scalar `"none"` disables all
transforms; `"identity"` is its compatibility alias. Original
measurements are never modified.

``` r

grouped <- mad_qc_bulk(bulk, metrics, transform = c(mapping_rate = "none"),
  group_by = "batch", verbose = FALSE)
attr(grouped, "mad_qc")$thresholds
#> # A tibble: 4 × 15
#>   batch metric   direction transform n_observations n_evaluated n_missing median
#>   <chr> <chr>    <chr>     <chr>              <int>       <int>     <int>  <dbl>
#> 1 a     library… lower     log1p                  3           3         0  17.2 
#> 2 b     library… lower     log1p                  2           2         0  17.3 
#> 3 a     mapping… lower     none                   3           3         0   0.92
#> 4 b     mapping… lower     none                   2           2         0   0.94
#> # ℹ 7 more variables: mad <dbl>, lower <dbl>, upper <dbl>, lower_raw <dbl>,
#> #   upper_raw <dbl>, n_outliers <int>, outlier_proportion <dbl>
```

## Inspect before filtering

``` r

flagged[flagged$library_size_mad_outlier %in% TRUE, ]
#>         library_size mapping_rate batch library_size_mad_outlier
#> sample1        1e+06          0.5     a                     TRUE
#>         mapping_rate_mad_outlier
#> sample1                     TRUE
```

This is inspection, not an automatic filtering recommendation. MAD
thresholds are adaptive heuristics rather than universal biological
cutoffs. Upper RNA count flags are not doublet calls, and doublet
detection is outside package scope. The playful console grade is
presentation only, never evidence of a biological conclusion or
laboratory failure. For scientific context, see the [scrapper
manual](https://bioconductor.org/packages/release/bioc/manuals/scrapper/man/scrapper.pdf).

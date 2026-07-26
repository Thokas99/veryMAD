# Summarize an existing MAD QC report

Summarizes flags already calculated by
[`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md).
Thresholds are not recalculated.

## Usage

``` r
summarize_mad_qc(qc, level = c("metric", "observation"), group_by = NULL)
```

## Arguments

- qc:

  A tidy report created by
  [`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md),
  either retaining its class or represented as an ordinary data frame
  with the required columns.

- level:

  Summarize by metric or original observation.

- group_by:

  `NULL` or one or more grouping column names in `qc`.

## Value

A tibble with deterministic columns and first-seen group, metric, or
observation order.

## Examples

``` r
d <- data.frame(sample = c("a", "a", "b"), x = c(1, 2, 10))
qc <- mad_qc(d, c(x = "upper"), group_by = "sample")
summarize_mad_qc(qc, "metric", group_by = "sample")
#> # A tibble: 2 × 8
#>   sample metric n_observations n_evaluated n_missing n_outliers
#>   <chr>  <chr>           <int>       <int>     <int>      <int>
#> 1 a      x                   2           2         0          0
#> 2 b      x                   1           1         0          0
#> # ℹ 2 more variables: outlier_proportion <dbl>, direction <chr>
summarize_mad_qc(qc, "observation")
#> # A tibble: 3 × 8
#>    .obs id    n_metrics n_evaluated n_missing n_outliers failed_metrics
#>   <int> <chr>     <int>       <int>     <int>      <int> <chr>         
#> 1     1 1             1           1         0          0 ""            
#> 2     2 2             1           1         0          0 ""            
#> 3     3 3             1           1         0          0 ""            
#> # ℹ 1 more variable: mad_qc_outlier <lgl>
```

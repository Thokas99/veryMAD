# Bulk observation-level MAD quality control

Rows are observations and `metrics` explicitly selects numeric QC
columns. A gene-by-sample expression matrix is not automatically a
QC-metric table. Selected metrics use `log1p` unless `transform` is
`"none"` or supplies a named partial override. Measurements are never
changed or filtered.

## Usage

``` r
mad_qc_bulk(
  data,
  metrics = NULL,
  nmads = 3,
  transform = "log1p",
  group_by = NULL,
  verbose = TRUE,
  overwrite = FALSE
)
```

## Arguments

- data:

  A matrix, data frame, or tibble with observations in rows.

- metrics:

  Named character vector mapping columns to `"lower"`, `"upper"`, or
  `"both"`. Required; veryMAD never guesses metrics.

- nmads:

  Positive number of MADs used for thresholds.

- transform:

  `"log1p"`, `"none"`, `"identity"`, or `"log10"`, or a named partial
  override. Unnamed metrics in an override use `"log1p"`.

- group_by:

  Optional unique column names for within-group calculations.

- verbose:

  Print a concise QC summary?

- overwrite:

  Replace existing per-metric flag columns?

## Value

The input as a data frame with one logical `*_mad_outlier` column per
metric. The compact report is available as `attr(result, "mad_qc")`.

## Examples

``` r
d <- data.frame(size = c(10, 11, 12, 100), rate = c(.9, .91, .92, .2))
out <- mad_qc_bulk(d, c(size = "upper", rate = "lower"),
                   transform = c(rate = "none"), verbose = FALSE)
attr(out, "mad_qc")$thresholds
#> # A tibble: 2 × 14
#>   metric direction transform n_observations n_evaluated n_missing median    mad
#>   <chr>  <chr>     <chr>              <int>       <int>     <int>  <dbl>  <dbl>
#> 1 size   upper     log1p                  4           4         0  2.52  0.124 
#> 2 rate   lower     none                   4           4         0  0.905 0.0148
#> # ℹ 6 more variables: lower <dbl>, upper <dbl>, lower_raw <dbl>,
#> #   upper_raw <dbl>, n_outliers <int>, outlier_proportion <dbl>
```

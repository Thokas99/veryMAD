# Single-cell observation-level MAD quality control

Accepts cell metadata or a Seurat object. Annotation adds only
`mad_qc_outlier`; individual flags remain in the compact report. `TRUE`
means at least one metric failed, `FALSE` means all passed, and `NA`
means none failed but at least one could not be evaluated. Cells are
never filtered.

## Usage

``` r
mad_qc_sc(
  object,
  metrics = NULL,
  nmads = 3,
  transform = "log1p",
  group_by = NULL,
  verbose = TRUE,
  overwrite = FALSE,
  action = c("annotate", "report")
)
```

## Arguments

- object:

  A data frame, tibble, or Seurat object.

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

- action:

  Annotate the input or return the compact report.

## Value

For `action = "report"`, a compact `mad_qc_result`. Otherwise the
metadata data frame or Seurat object with one `mad_qc_outlier` column.
The data-frame annotation stores the report in `attr(result, "mad_qc")`.

## Examples

``` r
cells <- data.frame(counts = c(NA, 10, 11, 100), mt = c(2, 3, 4, 30))
mad_qc_sc(cells, c(counts = "both", mt = "upper"),
          transform = c(mt = "none"), action = "report", verbose = FALSE)
#> $flags
#> # A tibble: 4 × 4
#>    .obs id    counts mt   
#>   <int> <chr> <lgl>  <lgl>
#> 1     1 1     NA     FALSE
#> 2     2 2     FALSE  FALSE
#> 3     3 3     FALSE  FALSE
#> 4     4 4     TRUE   TRUE 
#> 
#> $thresholds
#> # A tibble: 2 × 14
#>   metric direction transform n_observations n_evaluated n_missing median   mad
#>   <chr>  <chr>     <chr>              <int>       <int>     <int>  <dbl> <dbl>
#> 1 counts both      log1p                  4           3         1   2.48 0.129
#> 2 mt     upper     none                   4           4         0   3.5  1.48 
#> # ℹ 6 more variables: lower <dbl>, upper <dbl>, lower_raw <dbl>,
#> #   upper_raw <dbl>, n_outliers <int>, outlier_proportion <dbl>
#> 
#> $settings
#> $settings$metrics
#>  counts      mt 
#>  "both" "upper" 
#> 
#> $settings$transform
#>  counts      mt 
#> "log1p"  "none" 
#> 
#> $settings$nmads
#> [1] 3
#> 
#> $settings$group_by
#> NULL
#> 
#> $settings$observation_ids
#> [1] "1" "2" "3" "4"
#> 
#> 
#> attr(,"class")
#> [1] "mad_qc_result" "list"         
```

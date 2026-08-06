# Explicit MAD quality-control annotation

Rows are observations and `metrics` explicitly selects numeric QC
columns. One reference distribution is calculated per metric across all
observations; observations are never filtered automatically.

## Usage

``` r
mad_qc(
  data,
  metrics,
  nmads = 3,
  transform = "none",
  output = c("annotate", "report"),
  combine = TRUE,
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

  A non-empty named character vector mapping columns to `"lower"`,
  `"upper"`, or `"both"`.

- nmads:

  Positive number of MADs from the median for the limits.

- transform:

  A scalar transformation or named partial overrides. Supported values
  are `"none"`, `"log1p"`, and `"log10"`.

- output:

  Return annotated data (`"annotate"`) or a compact report (`"report"`).

- combine:

  Add the combined `mad_qc_outlier` flag?

- min_n:

  Minimum number of finite, non-missing observations required.

- zero_mad:

  Handling of a zero MAD: `"na"`, `"zero"`, or `"error"`.

- overwrite:

  Replace existing generated flag columns?

- verbose:

  Print a concise neutral summary?

## Value

Annotated input or a `verymad_qc` list with `flags`, `thresholds`, and
`settings`.

## Examples

``` r
x <- data.frame(low = c(1, 10, 11, 12, 13), high = c(1, 2, 3, 4, 20))
mad_qc(x, c(low = "lower", high = "upper"), verbose = FALSE)
#>   low high low_mad_outlier high_mad_outlier mad_qc_outlier
#> 1   1    1            TRUE            FALSE           TRUE
#> 2  10    2           FALSE            FALSE          FALSE
#> 3  11    3           FALSE            FALSE          FALSE
#> 4  12    4           FALSE            FALSE          FALSE
#> 5  13   20           FALSE             TRUE           TRUE
```

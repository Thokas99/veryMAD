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
qc_metadata <- data.frame(
  sample = paste0("sample_", 1:6),
  library_size = c(2.4e6, 2.5e6, 2.6e6, 2.7e6, 0.5e6, 2.5e6),
  pct_mito = c(.04, .05, .06, .05, .07, .30)
)
qc <- mad_qc(
  qc_metadata,
  metrics = c(library_size = "lower", pct_mito = "upper"),
  verbose = FALSE
)
qc[, c("sample", "mad_qc_outlier")]
#>     sample mad_qc_outlier
#> 1 sample_1          FALSE
#> 2 sample_2          FALSE
#> 3 sample_3          FALSE
#> 4 sample_4          FALSE
#> 5 sample_5           TRUE
#> 6 sample_6           TRUE
```

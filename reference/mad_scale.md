# Robust MAD scaling

Scale a numeric vector, matrix, or numeric data frame by its median and
median absolute deviation. For matrices, `margin = 1` scales rows and
`margin = 2` scales columns.

## Usage

``` r
mad_scale(
  x,
  center = TRUE,
  scale = TRUE,
  constant = 1.4826,
  na_rm = TRUE,
  zero_mad = c("zero", "na", "error"),
  margin = 2
)
```

## Arguments

- x:

  A numeric vector, matrix, or numeric data frame.

- center:

  Subtract the median?

- scale:

  Divide by the MAD?

- constant:

  Positive MAD consistency constant.

- na_rm:

  Remove missing values when calculating medians and MADs?

- zero_mad:

  Use zero, return `NA`, or error for a zero MAD.

- margin:

  Matrix margin, `1` for rows or `2` for columns. Names `"rows"` and
  `"columns"` are also accepted.

## Value

A numeric vector or matrix. Matrix input preserves dimensions and
dimnames; data-frame input returns a matrix.

## Examples

``` r
gene_expression <- c(gene_A = 8.2, gene_B = 8.7, gene_C = 18.0)
mad_scale(gene_expression)
#>     gene_A     gene_B     gene_C 
#> -0.6744908  0.0000000 12.5455281 

expression_matrix <- matrix(
  1:12,
  nrow = 3,
  dimnames = list(paste0("gene_", 1:3), paste0("sample_", 1:4))
)
mad_scale(expression_matrix, margin = 1)
#>         sample_1   sample_2  sample_3 sample_4
#> gene_1 -1.011736 -0.3372454 0.3372454 1.011736
#> gene_2 -1.011736 -0.3372454 0.3372454 1.011736
#> gene_3 -1.011736 -0.3372454 0.3372454 1.011736
```

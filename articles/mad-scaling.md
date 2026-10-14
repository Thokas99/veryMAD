# Robust MAD scaling

[`mad_scale()`](https://thokas99.github.io/veryMAD/reference/mad_scale.md)
centers numeric data by its median and divides by its median absolute
deviation (MAD). Set `center` or `scale` to `FALSE` to skip that
operation.

## Scale a vector

This example scales expression summaries for five genes.

``` r

gene_expression <- c(
  gene_A = 8.2,
  gene_B = 8.7,
  gene_C = 9.1,
  gene_D = 9.4,
  gene_E = 18.0
)

scaled_expression <- mad_scale(gene_expression)
round(scaled_expression, 2)
#> gene_A gene_B gene_C gene_D gene_E 
#>  -1.52  -0.67   0.00   0.51  15.01
```

Names and missing values are preserved.

## Scale an expression matrix

For a matrix, `margin = 1` scales rows and `margin = 2` scales columns.
Here rows are genes and columns are samples.

``` r

expression_matrix <- matrix(
  c(
    8.2, 8.5, 8.4, 8.6,
    9.1, 9.0, 9.3, 9.2,
    4.0, 4.2, NA, 4.1
  ),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(
    c("gene_A", "gene_B", "gene_C"),
    c("sample_1", "sample_2", "sample_3", "sample_4")
  )
)

scaled_by_sample <- mad_scale(expression_matrix, margin = 2)
round(scaled_by_sample, 2)
#>        sample_1 sample_2 sample_3 sample_4
#> gene_A     0.00     0.00    -0.67     0.00
#> gene_B     0.67     0.67     0.67     0.67
#> gene_C    -3.15    -5.80       NA    -5.06
```

Dimensions and dimnames are preserved. The default `zero_mad = "na"`
returns `NA` for zero-MAD groups. `"zero"` uses a unit denominator;
centered results are unscaled deviations. Choose `"error"` to stop on
zero MAD.

``` r

mostly_constant <- c(value_1 = 1, value_2 = 1, value_3 = 1, outlier = 100)
mad_scale(mostly_constant)
#> value_1 value_2 value_3 outlier 
#>      NA      NA      NA      NA
mad_scale(mostly_constant, zero_mad = "zero")
#> value_1 value_2 value_3 outlier 
#>       0       0       0      99
```

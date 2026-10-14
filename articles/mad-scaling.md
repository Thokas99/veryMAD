# Robust MAD scaling

[`mad_scale()`](https://thokas99.github.io/veryMAD/reference/mad_scale.md)
optionally subtracts the median and divides by the scaled median
absolute deviation (MAD). The default MAD constant, 1.4826, is
consistent with the standard deviation for normally distributed data.

## Vectors

``` r

gene_expression <- c(
  gene_A = 8.2,
  gene_B = 8.7,
  gene_C = 9.1,
  gene_D = 9.4,
  gene_E = 18.0
)

round(mad_scale(gene_expression), 2)
#> gene_A gene_B gene_C gene_D gene_E 
#>  -1.52  -0.67   0.00   0.51  15.01
```

Names and missing values are preserved.

## Matrices

For matrices, the margin selects which groups are scaled:

| `margin`           | Groups scaled |
|--------------------|---------------|
| `1` or `"rows"`    | Each row      |
| `2` or `"columns"` | Each column   |

Here, rows are genes and columns are samples.

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

round(mad_scale(expression_matrix, margin = "columns"), 2)
#>        sample_1 sample_2 sample_3 sample_4
#> gene_A     0.00     0.00    -0.67     0.00
#> gene_B     0.67     0.67     0.67     0.67
#> gene_C    -3.15    -5.80       NA    -5.06
```

Dimensions and dimnames are preserved. When a group has zero MAD, choose
how to handle it:

| `zero_mad` | Behavior when `scale = TRUE` |
|----|----|
| `"na"` (default) | Return `NA` for the affected group |
| `"zero"` | Use a unit denominator; centered values are deviations from the median |
| `"error"` | Stop |

``` r

constant_values <- c(value_1 = 1, value_2 = 1, value_3 = 1, outlier = 100)
mad_scale(constant_values)
#> value_1 value_2 value_3 outlier 
#>      NA      NA      NA      NA
mad_scale(constant_values, zero_mad = "zero")
#> value_1 value_2 value_3 outlier 
#>       0       0       0      99
```

For QC flags and their visual review, see [Getting
started](https://thokas99.github.io/veryMAD/articles/getting-started.qmd)
and [QC
plots](https://thokas99.github.io/veryMAD/articles/qc-plots.qmd).

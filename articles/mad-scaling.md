# Robust MAD scaling

[`mad_scale()`](https://thokas99.github.io/veryMAD/reference/mad_scale.md)
centers numeric data by its median and, by default, divides by the
median absolute deviation (MAD). This is useful when a few extreme
values should not dominate the scale used for an exploratory analysis.

## Scale a vector

This small vector represents expression summaries for five genes. The
largest value is easy to see, but it should not determine the centre and
scale by itself.

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

The names are preserved, and missing values remain missing. Set
`center = FALSE` or `scale = FALSE` when only one part of the
transformation is needed.

## Scale an expression matrix

For a matrix, `margin = 1` scales rows and `margin = 2` scales columns.
Here, rows are genes and columns are samples, so column scaling puts
samples on a common robust scale.

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

The matrix dimensions and dimnames are preserved. A zero MAD can be
handled with `zero_mad = "zero"`, `"na"`, or `"error"`, depending on how
the result should be used downstream.

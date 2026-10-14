# Center and scale with the median and MAD

Optionally center on the median and scale by the median absolute
deviation (MAD). For matrices, `margin = 1` operates on rows and
`margin = 2` on columns.

## Usage

``` r
mad_scale(
  x,
  center = TRUE,
  scale = TRUE,
  constant = 1.4826,
  na_rm = TRUE,
  zero_mad = c("na", "zero", "error"),
  margin = 2
)
```

## Arguments

- x:

  Numeric vector, matrix, or data frame with numeric columns.

- center:

  Subtract the median from each vector or matrix margin?

- scale:

  Divide each vector or matrix margin by its MAD?

- constant:

  Positive consistency constant applied to the MAD. The default, 1.4826,
  makes it consistent with the standard deviation under a normal
  distribution.

- na_rm:

  Omit missing values when estimating medians and MADs. Missing input
  values remain missing in the result.

- zero_mad:

  Policy for a zero MAD when `scale = TRUE`: `"na"` returns `NA` for the
  affected vector or margin, `"zero"` uses a unit denominator, and
  `"error"` stops. Centered results under `"zero"` are deviations from
  the median without MAD scaling.

- margin:

  Matrix margin: `1` or `"rows"` operates on rows; `2` or `"columns"`
  operates on columns. Ignored for vectors.

## Value

A numeric vector or matrix. Matrices preserve dimensions and dimnames;
data frames are returned as matrices.

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

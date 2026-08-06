# Robust MAD scaling

Scale a numeric vector, matrix, or numeric data frame by its median and
median absolute deviation. For matrices, `margin = 1` scales rows and
`margin = 2` scales columns. `mad_z_score()` is an alias.

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

mad_z_score(...)
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

- ...:

  Additional arguments are passed to `mad_scale()`.

## Value

A numeric vector or matrix. Matrix input preserves dimensions and
dimnames; data-frame input returns a matrix.

## Examples

``` r
mad_scale(c(a = 1, b = 2, c = 100))
#>          a          b          c 
#> -0.6744908  0.0000000 66.1000944 
mad_scale(matrix(1:12, nrow = 3), margin = 1)
#>           [,1]       [,2]      [,3]     [,4]
#> [1,] -1.011736 -0.3372454 0.3372454 1.011736
#> [2,] -1.011736 -0.3372454 0.3372454 1.011736
#> [3,] -1.011736 -0.3372454 0.3372454 1.011736
```

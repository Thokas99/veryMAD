# Robust MAD scaling for vectors and dense matrices

`mad_scale()` median-centres and MAD-scales a numeric vector or each
margin of an ordinary dense numeric matrix. For expression heatmaps, use
`margin = "rows"` and pass the result to the plotting function without
applying another scaling step.

## Usage

``` r
mad_scale(
  x,
  margin = c("rows", "columns"),
  constant = 1.4826,
  na_rm = TRUE,
  zero_mad = c("zero", "na", "error")
)
```

## Arguments

- x:

  A numeric vector or ordinary dense numeric matrix.

- margin:

  Matrix margin to scale. Ignored for vector input.

- constant:

  Positive MAD consistency constant.

- na_rm:

  Remove missing values when calculating centres and spreads?

- zero_mad:

  Behavior for zero-MAD vectors or matrix margins.

## Value

A numeric vector or matrix preserving names or dimnames.

## Examples

``` r
mad_scale(c(a = 1, b = 2, c = 100))
#>          a          b          c 
#> -0.6744908  0.0000000 66.1000944 
x <- matrix(1:12, nrow = 3, dimnames = list(paste0("gene", 1:3), NULL))
mad_scale(x, margin = "rows")
#>            [,1]       [,2]      [,3]     [,4]
#> gene1 -1.011736 -0.3372454 0.3372454 1.011736
#> gene2 -1.011736 -0.3372454 0.3372454 1.011736
#> gene3 -1.011736 -0.3372454 0.3372454 1.011736
```

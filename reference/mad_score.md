# Calculate robust MAD scores

Calculate robust MAD scores

## Usage

``` r
mad_score(
  x,
  constant = 1.4826,
  na_rm = TRUE,
  zero_mad = c("zero", "na", "error")
)
```

## Arguments

- x:

  A numeric vector without infinite values.

- constant:

  Positive MAD consistency constant.

- na_rm:

  Remove missing values before calculation?

- zero_mad:

  Behavior when MAD is zero: neutral limits, missing limits, or an
  error.

## Value

A numeric vector with the length and names of `x`.

## Examples

``` r
mad_score(c(a = 1, b = 2, c = 100))
#>          a          b          c 
#> -0.6744908  0.0000000 66.1000944 
```

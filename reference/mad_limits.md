# Calculate MAD-based limits

Calculate MAD-based limits

## Usage

``` r
mad_limits(
  x,
  nmads = 3,
  direction = c("both", "lower", "upper"),
  constant = 1.4826,
  na_rm = TRUE,
  zero_mad = c("zero", "na", "error")
)
```

## Arguments

- x:

  A numeric vector without infinite values.

- nmads:

  Number of MADs from the median.

- direction:

  Tail(s) for which to calculate finite limits.

- constant:

  Positive MAD consistency constant.

- na_rm:

  Remove missing values before calculation?

- zero_mad:

  Behavior when MAD is zero: neutral limits, missing limits, or an
  error.

## Value

A one-row tibble with `median`, `mad`, `lower`, and `upper`.

## Examples

``` r
mad_limits(c(1, 2, 2, 3, 100), direction = "upper")
#> # A tibble: 1 × 4
#>   median   mad lower upper
#>    <dbl> <dbl> <dbl> <dbl>
#> 1      2  1.48    NA  6.45
```

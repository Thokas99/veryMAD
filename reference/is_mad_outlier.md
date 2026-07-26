# Flag MAD outliers

Flag MAD outliers

## Usage

``` r
is_mad_outlier(
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

A logical vector with the length and names of `x`.

## Examples

``` r
is_mad_outlier(c(1, 2, 2, 3, 100), direction = "upper")
#> [1] FALSE FALSE FALSE FALSE  TRUE
```

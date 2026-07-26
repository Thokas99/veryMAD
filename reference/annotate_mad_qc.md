# Add MAD QC flags to a data frame

Add MAD QC flags to a data frame

## Usage

``` r
annotate_mad_qc(
  data,
  qc,
  suffix = "_mad_outlier",
  overall = "mad_qc_outlier",
  overwrite = FALSE
)
```

## Arguments

- data:

  Original data frame used by
  [`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md).

- qc:

  A tidy result from
  [`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md).

- suffix:

  Suffix for per-metric logical columns.

- overall:

  Name of the overall logical flag.

- overwrite:

  Replace only conflicting veryMAD flag columns?

## Value

`data` with aligned logical QC columns appended.

## Examples

``` r
d <- data.frame(x = c(1, 1, 1, 10))
q <- mad_qc(d, c(x = "upper"))
annotate_mad_qc(d, q)
#>    x x_mad_outlier mad_qc_outlier
#> 1  1         FALSE          FALSE
#> 2  1         FALSE          FALSE
#> 3  1         FALSE          FALSE
#> 4 10         FALSE          FALSE
```

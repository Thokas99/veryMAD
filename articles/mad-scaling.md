# Robust MAD scaling

``` r

library(veryMAD)
mad_scale(c(1, 2, 100))
```

    ## [1] -0.6744908  0.0000000 66.1000944

For matrices, `margin = 1` scales rows and `margin = 2` scales columns.
Matrix dimensions and names are preserved; numeric data frames are
returned as matrices. A MAD-based scale is robust to extreme values and
differs from a standard deviation-based z-score. Zero-MAD margins can
return zeros, `NA`, or an error through `zero_mad`.

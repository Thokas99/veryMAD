# Getting started with veryMAD

`veryMAD` accepts observation metadata and requires users to select
numeric metrics and directions explicitly.

``` r

metadata <- data.frame(
  library_size = c(2e6, 31e6, 32e6, 30e6, 33e6),
  mapping_rate = c(.50, .92, .94, .93, .95),
  row.names = paste0("sample", 1:5)
)
mad_qc(metadata, c(library_size = "lower", mapping_rate = "lower"),
  transform = c(library_size = "log1p"), verbose = FALSE)
#>         library_size mapping_rate library_size_mad_outlier
#> sample1      2.0e+06         0.50                     TRUE
#> sample2      3.1e+07         0.92                    FALSE
#> sample3      3.2e+07         0.94                    FALSE
#> sample4      3.0e+07         0.93                    FALSE
#> sample5      3.3e+07         0.95                    FALSE
#>         mapping_rate_mad_outlier mad_qc_outlier
#> sample1                     TRUE           TRUE
#> sample2                    FALSE          FALSE
#> sample3                    FALSE          FALSE
#> sample4                    FALSE          FALSE
#> sample5                    FALSE          FALSE
```

Use `output = "report"` to inspect one threshold row per metric and one
flag row per observation. Missing values remain `NA`; no observations
are removed.

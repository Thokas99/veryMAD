# Robust MAD scaling

``` r

set.seed(303)
gene_scores <- rlnorm(120, meanlog = 4, sdlog = .8)
names(gene_scores) <- sprintf("gene_%03d", seq_along(gene_scores))
scaled_scores <- mad_scale(gene_scores)
head(round(scaled_scores, 2))
#> gene_001 gene_002 gene_003 gene_004 gene_005 gene_006 
#>    -0.41     0.29     0.42     0.06    -0.46     2.25
```

The vector example has 120 values, allowing the robust centre and scale
to be estimated from a more realistic feature set. For matrices,
`margin = 1` scales rows and `margin = 2` scales columns. The next
example uses a 200 by 24 expression matrix and includes one missing
value to demonstrate that missing values remain missing.

``` r

expression <- matrix(
  rlnorm(200 * 24, meanlog = 4, sdlog = .8),
  nrow = 200,
  ncol = 24,
  dimnames = list(
    sprintf("gene_%03d", 1:200),
    sprintf("sample_%02d", 1:24)
  )
)
expression[, c(4, 19)] <- expression[, c(4, 19)] * c(2.5, .5)
expression[12, 7] <- NA_real_

scaled_genes <- mad_scale(expression, margin = 1)
scaled_samples <- mad_scale(expression, margin = 2)

dim(scaled_genes)
#> [1] 200  24
dimnames(scaled_genes)[[1]][1:3]
#> [1] "gene_001" "gene_002" "gene_003"
round(scaled_genes[1:5, 1:5], 2)
#>          sample_01 sample_02 sample_03 sample_04 sample_05
#> gene_001      3.55     -0.74      0.20      6.43      1.14
#> gene_002      0.50     -1.42      0.01     -0.73      0.60
#> gene_003     -0.39     -0.19     -0.39      1.88      4.14
#> gene_004      1.63     -0.18      0.72      1.84     -0.26
#> gene_005     -0.63      1.38      1.75      0.26     -0.76
round(scaled_samples[1:5, 1:5], 2)
#>          sample_01 sample_02 sample_03 sample_04 sample_05
#> gene_001      3.14     -0.81     -0.13      2.97      0.88
#> gene_002      0.45     -1.03     -0.12     -0.45      0.48
#> gene_003     -0.41     -0.22     -0.58      1.13      5.19
#> gene_004      0.96     -0.58      0.03      0.44     -0.66
#> gene_005     -0.49      2.16      2.98      0.29     -0.72
is.na(scaled_genes[12, 7])
#> [1] TRUE
```

Dimensions and dimnames are preserved; numeric data frames are returned
as matrices. A MAD-based scale is robust to extreme values and differs
from a standard-deviation-based scaling. Zero-MAD margins can return
zeros, `NA`, or an error through `zero_mad`.

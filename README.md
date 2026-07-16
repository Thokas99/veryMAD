# veryMAD

Robust median absolute deviation (MAD) helpers for QC metrics and
expression-like matrices in R.

veryMAD is useful when a few extreme observations should not control the
center, scale, or feature ranking. It provides small functions for MAD scores,
outlier flags, winsorization, robust matrix scaling, QC flagging, and variable
feature selection.

## Installation

Install from a local checkout:

```sh
R CMD INSTALL .
```

Or install a built source package:

```r
install.packages("veryMAD_0.0.0.9000.tar.gz", repos = NULL, type = "source")
```

## Quick Start

```r
library(veryMAD)

x <- c(1, 2, 2, 3, 100)

mad_score(x)
is_outlier(x)
winsorize(x)
```

`mad_score()` returns robust z-like scores using the median and MAD.
`is_outlier()` flags values whose absolute MAD score is above the threshold
(`3.5` by default). `winsorize()` caps those values at the MAD-based limits.

## Matrix Helpers

```r
mat <- matrix(
  c(1, 2, 2, 3, 100,
    10, 20, 20, 30, 40),
  nrow = 2,
  byrow = TRUE
)
rownames(mat) <- c("gene_a", "gene_b")

row_mad(mat)
robust_scale(mat)
select_variable_features(mat, n = 1)
```

Use `margin = 1` for rows and `margin = 2` for columns.

## QC Tables

```r
qc <- data.frame(
  nCount_RNA = c(500, 600, 550, 20000, 580),
  percent_mt = c(2, 3, 2.5, 3, 40)
)

mad_qc(qc, metrics = c("nCount_RNA", "percent_mt"))
```

`mad_qc()` appends one logical outlier column per metric plus a combined
`mad_qc_outlier` flag.

## API

| Function | Purpose |
| --- | --- |
| `mad_score()` | Robust z-like score from median and MAD. |
| `is_outlier()` | Logical MAD-threshold outlier flag. |
| `winsorize()` | Cap values at MAD-threshold limits. |
| `row_mad()` | Row-wise or column-wise MAD for numeric matrices. |
| `robust_scale()` | Median-center and MAD-scale numeric matrices. |
| `mad_qc()` | Add MAD-based QC flags to data frames or Seurat metadata. |
| `select_variable_features()` | Select rows or columns with highest MAD. |

## Development

```r
R CMD build .
R CMD check --no-manual --no-build-vignettes veryMAD_0.0.0.9000.tar.gz
```

# veryMAD

Robust median absolute deviation (MAD) helpers for QC metrics and
expression-like matrices in R.

veryMAD is useful when a few extreme observations should not control the
center, scale, or feature ranking. It provides small functions for MAD scores,
outlier flags, winsorization, robust matrix scaling, QC flagging, and variable
feature selection.

## Installation

Install the development version from GitHub with pak:

```r
install.packages("pak")
pak::pak("Thokas99/veryMAD")
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

mad_qc(qc, metrics = c(nCount_RNA = "both", percent_mt = "upper"))
```

`mad_qc()` expects a named character vector: names are metric columns and
values are outlier directions (`"lower"`, `"upper"`, or `"both"`).

It returns a tidy long data frame with one row per observation and metric:
`id`, `metric`, `value`, `median`, `mad`, `lower`, `upper`, and `is_outlier`.

## Verbose Output

All exported functions are quiet by default. Set `verbose = TRUE` to get a
short `cli` summary without changing the returned object.

```r
is_outlier(x, verbose = TRUE)
mad_qc(qc, metrics = c(nCount_RNA = "both", percent_mt = "upper"), verbose = TRUE)
```

`mad_qc(verbose = TRUE)` also shows a progress bar while it checks metrics.

## API

| Function | Purpose |
| --- | --- |
| `mad_score()` | Robust z-like score from median and MAD. |
| `is_outlier()` | Logical MAD-threshold outlier flag. |
| `winsorize()` | Cap values at MAD-threshold limits. |
| `row_mad()` | Row-wise or column-wise MAD for numeric matrices. |
| `robust_scale()` | Median-center and MAD-scale numeric matrices. |
| `mad_qc()` | Return tidy MAD-based QC thresholds and outlier flags. |
| `select_variable_features()` | Select rows or columns with highest MAD. |

## Changelog

Project changes are tracked in [CHANGELOG.md](CHANGELOG.md). The changelog
format follows [Keep a Changelog 1.1.0](https://keepachangelog.com/en/1.1.0/).

## Development

```r
R CMD build .
R CMD check --no-manual --no-build-vignettes veryMAD_*.tar.gz
```

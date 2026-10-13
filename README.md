<p align="center">
  <img src="man/figures/veryMAD-logo.svg" alt="veryMAD logo" width="240">
</p>

<h1 align="center">veryMAD</h1>

<p align="center"><strong>MAD quality control and robust scaling for R</strong></p>

<p align="center">
  <a href="https://github.com/Thokas99/veryMAD/actions/workflows/R-CMD-check.yaml"><img src="https://github.com/Thokas99/veryMAD/actions/workflows/R-CMD-check.yaml/badge.svg" alt="R CMD check"></a>
  <a href="https://github.com/Thokas99/veryMAD/actions/workflows/pkgdown.yaml"><img src="https://github.com/Thokas99/veryMAD/actions/workflows/pkgdown.yaml/badge.svg" alt="pkgdown"></a>
  <a href="https://github.com/Thokas99/veryMAD/releases"><img src="https://img.shields.io/github/v/release/Thokas99/veryMAD?display_name=tag&sort=semver" alt="Latest release"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-blue.svg" alt="MIT license"></a>
</p>

veryMAD computes median absolute deviation (MAD) thresholds for selected
numeric metrics. It adds observation-level flags or returns a threshold report.
Rows stay in the output; filtering remains part of your analysis.

## Install

```r
pak::pak("Thokas99/veryMAD")
```

## Annotate QC metrics

Choose each metric, its outlier direction, and any transformation. For example,
count-like RNA metrics can use an explicit `log1p` transformation.

```r
library(veryMAD)

metadata <- data.frame(
  nCount_RNA = c(820, 910, 970, 1050, 1120, 100),
  nFeature_RNA = c(430, 480, 520, 570, 610, 55),
  percent.mt = c(5, 6, 7, 8, 9, 35)
)

annotated <- mad_qc(
  metadata,
  metrics = c(nCount_RNA = "lower", nFeature_RNA = "lower",
              percent.mt = "upper"),
  transform = c(nCount_RNA = "log1p", nFeature_RNA = "log1p"),
  verbose = FALSE
)

annotated[, c("nCount_RNA_mad_outlier", "nFeature_RNA_mad_outlier",
              "percent.mt_mad_outlier", "mad_qc_outlier")]
```

For Seurat metadata, `percent.mt` is usually stored on a 0–100 scale: `8`
means 8%.

## Review thresholds

Set `output = "report"` to return flags, thresholds, and settings without
adding columns to the input.

```r
report <- mad_qc(
  metadata,
  metrics = c(nCount_RNA = "lower", nFeature_RNA = "lower",
              percent.mt = "upper"),
  transform = c(nCount_RNA = "log1p", nFeature_RNA = "log1p"),
  output = "report",
  verbose = FALSE
)

report$thresholds[, c("metric", "direction", "lower_raw", "upper_raw",
                      "status", "n_used", "n_flagged", "n_unresolved")]
```

The report returns flags, thresholds, and settings. Status identifies missing,
insufficient, and zero-MAD metrics; counts show evaluated and unresolved
observations. Zero-MAD flags are `NA` by default.

## Scale with MAD

```r
gene_expression <- c(gene_A = 8.2, gene_B = 8.7, gene_C = 18.0)
mad_scale(gene_expression)
```

For matrices, `margin = 1` scales rows and `margin = 2` scales columns.
Dimensions and dimnames are preserved. The default `zero_mad = "na"` returns
`NA` for zero-MAD groups. `zero_mad = "zero"` uses a unit denominator and
returns unscaled deviations when centered; `zero_mad = "error"` stops.

## Guides

- [Getting started](https://thokas99.github.io/veryMAD/articles/getting-started.html)
- [Interpreting QC](https://thokas99.github.io/veryMAD/articles/qc-interpretation.html)
- [MAD scaling](https://thokas99.github.io/veryMAD/articles/mad-scaling.html)
- [Function reference](https://thokas99.github.io/veryMAD/reference/)
- [Release notes](NEWS.md)

## Further reading

- [Subramanian et al. (2022), adaptive QC in single-cell transcriptomics](https://doi.org/10.1186/s13059-022-02820-w)
- [Single-cell best practices: quality control](https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html)
- [OSCA: quality control](https://bioconductor.org/books/3.13/OSCA.basic/quality-control.html#common-choices-of-qc-metrics)
- Documentation guidance: [Kassis et al. (2026), *Scientific Agent Skills*](https://doi.org/10.48550/arXiv.2609.00065)

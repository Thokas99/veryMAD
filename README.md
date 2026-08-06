<p align="center">
  <img src="man/figures/veryMAD-logo.svg" alt="veryMAD logo" width="280">
</p>

<h1 align="center">veryMAD</h1>

<p align="center"><strong>Explicit MAD quality control</strong></p>

<p align="center">
  <a href="https://github.com/Thokas99/veryMAD/actions/workflows/R-CMD-check.yaml"><img src="https://github.com/Thokas99/veryMAD/actions/workflows/R-CMD-check.yaml/badge.svg" alt="R CMD check"></a>
  <a href="https://github.com/Thokas99/veryMAD/actions/workflows/pkgdown.yaml"><img src="https://github.com/Thokas99/veryMAD/actions/workflows/pkgdown.yaml/badge.svg" alt="pkgdown"></a>
  <a href="https://github.com/Thokas99/veryMAD/releases"><img src="https://img.shields.io/github/v/release/Thokas99/veryMAD?display_name=tag&sort=semver" alt="Latest release"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-blue.svg" alt="MIT license"></a>
</p>

> [!NOTE]
> veryMAD calculates transparent MAD thresholds and flags. It never filters
> observations, guesses metrics, or turns a statistical flag into a biological
> diagnosis.

## What veryMAD does

| Capability | Function | Output |
| --- | --- | --- |
| Observation-level QC | `mad_qc()` | Annotated metadata or a compact report |
| Robust scaling | `mad_scale()` | MAD-scaled vectors and matrices |

Rows are observations. Columns are QC metrics selected by you. Directions and
transformations are explicit, deterministic, and auditable.

## Install

Install veryMAD from GitHub with [`pak`](https://pak.r-lib.org/):

```r
pak::pak("Thokas99/veryMAD")
```

## MAD QC in one call

```r
library(veryMAD)

metadata <- data.frame(
  library_size = c(2e6, 31e6, 32e6, 30e6, 33e6),
  mapping_rate = c(0.50, 0.92, 0.94, 0.93, 0.95),
  row.names = paste0("sample", 1:5)
)

annotated <- mad_qc(
  metadata,
  metrics = c(
    library_size = "lower",
    mapping_rate = "lower"
  ),
  transform = c(library_size = "log1p")
)

annotated
```

The default annotation keeps the original data and adds one logical flag per
metric plus `mad_qc_outlier`:

```text
library_size_mad_outlier
mapping_rate_mad_outlier
mad_qc_outlier
```

> [!TIP]
> Use `transform = "none"` for raw-scale calculations, or provide named
> partial overrides such as `c(library_size = "log1p")`.

## Inspect thresholds with a report

```r
report <- mad_qc(
  metadata,
  metrics = c(library_size = "lower", mapping_rate = "lower"),
  transform = c(library_size = "log1p"),
  output = "report",
  verbose = FALSE
)

report$thresholds
report$flags
report$settings
```

The report contains only three components:

```r
names(report)
#> "flags" "thresholds" "settings"
```

Thresholds retain both calculation-scale and raw-scale limits. Missing,
undersized, and zero-MAD calculations remain visible through `status` rather
than being silently converted into a decision.

## Robust MAD scaling

```r
values <- c(a = 1, b = 2, c = 100)
mad_scale(values)

matrix_data <- matrix(1:12, nrow = 3,
  dimnames = list(paste0("gene", 1:3), paste0("sample", 1:4)))

mad_scale(matrix_data, margin = 1) # scale rows
mad_scale(matrix_data, margin = 2) # scale columns
```

`margin = 1` scales rows; `margin = 2` scales columns. Matrix dimensions and
dimnames are preserved. Numeric data frames are returned as matrices.

## Interpretation and guardrails

> [!WARNING]
> MAD thresholds are adaptive statistical heuristics. A flag is not a sample
> rejection, a laboratory diagnosis, or a doublet call.

- Metrics, directions, and transformations are always selected explicitly.
- `lower`, `upper`, and `both` describe statistical tails only.
- `min_n` is a computational safeguard, not a biological rule.
- Missing values stay missing.
- veryMAD does not filter data automatically.
- For stratified QC, split the data yourself and call `mad_qc()` separately.

## Documentation

- **[Package website](https://thokas99.github.io/veryMAD/)**
- **[Function reference](https://thokas99.github.io/veryMAD/reference/)**
- **[Getting started](https://thokas99.github.io/veryMAD/articles/getting-started.html)**
- **[QC interpretation](https://thokas99.github.io/veryMAD/articles/qc-interpretation.html)**
- **[MAD scaling guide](https://thokas99.github.io/veryMAD/articles/mad-scaling.html)**
- **[Release notes](NEWS.md)**

## Bioconductor and single-cell QC guidance

These sources motivate veryMAD's explicit MAD-based QC scope while also
supporting its guardrails: MAD flags are adaptive statistical evidence, not
automatic biological decisions. Thresholds should be reviewed alongside
multiple QC metrics, diagnostic plots, and the biological context of the
dataset.

- **Subramanian, A., Alperovich, M., Yang, Y., & Li, B. (2022).** Biology-inspired data-driven quality control for scientific discovery in single-cell transcriptomics. *Genome Biology, 23*, 267. https://doi.org/10.1186/s13059-022-02820-w — Presents adaptive MAD-based QC across cell-type or cluster contexts and emphasizes retaining biologically meaningful populations that fixed, data-agnostic thresholds can remove. [Article](https://pmc.ncbi.nlm.nih.gov/articles/PMC9793662/)
- **Single-Cell Best Practices Consortium. (n.d.).** Quality control. In *Single-cell best practices*. Retrieved August 6, 2026, from https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html — Recommends lenient MAD-based cutoffs for QC while warning that aggressive filtering can bias against smaller subpopulations and that several covariates should be considered together.
- **Amezquita, R. A., Lun, A. T. L., Hicks, S. C., & Gottardo, R. (2021).** Quality control. In *Basics of single-cell analysis with Bioconductor* (Bioconductor 3.13). Bioconductor. https://bioconductor.org/books/3.13/OSCA.basic/quality-control.html#common-choices-of-qc-metrics — Defines common QC metrics such as library size, detected features, and mitochondrial or spike-in proportions, and describes adaptive MAD thresholds with explicit direction and transformation choices.

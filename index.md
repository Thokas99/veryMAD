# veryMAD

![veryMAD logo](reference/figures/veryMAD-logo.svg)

**Explicit MAD quality control · robust scaling · no hidden decisions**

[![R CMD
check](https://github.com/Thokas99/veryMAD/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Thokas99/veryMAD/actions/workflows/R-CMD-check.yaml)
[![pkgdown](https://github.com/Thokas99/veryMAD/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/Thokas99/veryMAD/actions/workflows/pkgdown.yaml)
[![Latest
release](https://img.shields.io/github/v/release/Thokas99/veryMAD?display_name=tag&sort=semver)](https://github.com/Thokas99/veryMAD/releases)
[![MIT
license](https://img.shields.io/badge/license-MIT-blue.svg)](https://thokas99.github.io/veryMAD/LICENSE)

> \[!NOTE\] veryMAD calculates transparent MAD thresholds and flags. It
> never filters observations, guesses metrics, or turns a statistical
> flag into a biological diagnosis.

## What veryMAD does

| Capability | Function | Output |
|----|----|----|
| Observation-level QC | [`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md) | Annotated metadata or a compact report |
| Robust scaling | [`mad_scale()`](https://thokas99.github.io/veryMAD/reference/mad_scale.md) | MAD-scaled vectors and matrices |
| Familiar alias | [`mad_z_score()`](https://thokas99.github.io/veryMAD/reference/mad_scale.md) | Alias of [`mad_scale()`](https://thokas99.github.io/veryMAD/reference/mad_scale.md) |

Rows are observations. Columns are QC metrics selected by you.
Directions and transformations are explicit, deterministic, and
auditable.

## Install

Install the released version from GitHub with
[`pak`](https://pak.r-lib.org/):

``` r

pak::pak("Thokas99/veryMAD")
```

Or install the development version:

``` r

pak::pak("Thokas99/veryMAD@main")
```

## MAD QC in one call

``` r

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

The default annotation keeps the original data and adds one logical flag
per metric plus `mad_qc_outlier`:

``` text
library_size_mad_outlier
mapping_rate_mad_outlier
mad_qc_outlier
```

> \[!TIP\] Use `transform = "none"` for raw-scale calculations, or
> provide named partial overrides such as `c(library_size = "log1p")`.
> veryMAD never infers a transformation from a column name.

## Inspect thresholds with a report

``` r

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

``` r

names(report)
#> "flags" "thresholds" "settings"
```

Thresholds retain both calculation-scale and raw-scale limits. Missing,
undersized, and zero-MAD calculations remain visible through `status`
rather than being silently converted into a decision.

## Robust MAD scaling

``` r

values <- c(a = 1, b = 2, c = 100)
mad_scale(values)

matrix_data <- matrix(1:12, nrow = 3,
  dimnames = list(paste0("gene", 1:3), paste0("sample", 1:4)))

mad_scale(matrix_data, margin = 1) # scale rows
mad_scale(matrix_data, margin = 2) # scale columns
mad_z_score(matrix_data, margin = 1) # same implementation
```

`margin = 1` scales rows; `margin = 2` scales columns. Matrix dimensions
and dimnames are preserved. Numeric data frames are returned as
matrices.

## Interpretation and guardrails

> \[!WARNING\] MAD thresholds are adaptive statistical heuristics. A
> flag is not a sample rejection, a laboratory diagnosis, or a doublet
> call.

- Metrics, directions, and transformations are always selected
  explicitly.
- `lower`, `upper`, and `both` describe statistical tails only.
- `min_n` is a computational safeguard, not a biological rule.
- Missing values stay missing.
- veryMAD does not filter data automatically.
- For stratified QC, split the data yourself and call
  [`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md)
  separately.

## Documentation

- **[Package website](https://thokas99.github.io/veryMAD/)**
- **[Function
  reference](https://thokas99.github.io/veryMAD/reference/)**
- **[Getting
  started](https://thokas99.github.io/veryMAD/articles/getting-started.html)**
- **[QC
  interpretation](https://thokas99.github.io/veryMAD/articles/qc-interpretation.html)**
- **[MAD scaling
  guide](https://thokas99.github.io/veryMAD/articles/mad-scaling.html)**
- **[Release notes](https://thokas99.github.io/veryMAD/NEWS.md)**

## Scope of the 0.5.0 release

Version 0.5.0 is a deliberate breaking simplification. The public API is
now small enough to audit:
[`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md) for
observation-level QC and
[`mad_scale()`](https://thokas99.github.io/veryMAD/reference/mad_scale.md)
for robust scaling. Legacy grouped, modality-specific, plotting,
grading, and low-level statistical entry points are no longer public
functions.

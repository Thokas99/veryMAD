# Changelog

## [veryMAD (development version)](https://github.com/Thokas99/veryMAD/compare/v0.7.0...HEAD)

### Added

- Added raw-scale MAD cutoff overlays to
  [`plot_mad_qc_scatter()`](https://thokas99.github.io/veryMAD/reference/plot_mad_qc_scatter.md)
  when given a matching
  [`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md)
  report.
- Added a 10,000-cell SPARSim vignette with a composed, batch-aware QC
  plot.

### Changed

- Reduced point size and opacity to keep dense QC plots readable.

## [veryMAD 0.7.0](https://github.com/Thokas99/veryMAD/compare/v0.6.0...v0.7.0) - 2026-10-14

### Added

- Added batch-specific and reference-set thresholds to
  [`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md).
- Added
  [`plot_mad_qc()`](https://thokas99.github.io/veryMAD/reference/plot_mad_qc.md)
  for batch distributions and
  [`plot_mad_qc_scatter()`](https://thokas99.github.io/veryMAD/reference/plot_mad_qc_scatter.md)
  for paired metrics with per-batch Spearman rho and finite-pair counts.

### Changed

- Clarified Seurat metric definitions and QC review guidance.

## [veryMAD 0.6.0](https://github.com/Thokas99/veryMAD/compare/v0.5.1...v0.6.0) - 2026-10-14

### Added

- Added threshold counts and report-only output to
  [`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md).
- Added single-cell examples and a reproducible benchmark.

### Changed

- Set
  [`mad_scale()`](https://thokas99.github.io/veryMAD/reference/mad_scale.md)
  to `zero_mad = "na"` by default; `"zero"` uses a unit denominator and
  returns unscaled deviations when centered.
- Reserved `id` and `mad_qc_outlier` as metric names in
  [`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md).

## [veryMAD 0.5.1](https://github.com/Thokas99/veryMAD/compare/v0.5.0...v0.5.1) - 2026-08-06

### Added

- Established
  [`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md)
  and
  [`mad_scale()`](https://thokas99.github.io/veryMAD/reference/mad_scale.md)
  as the public API.
- Added the combined `mad_qc_outlier` flag and compact report printing.

## [veryMAD 0.5.0](https://github.com/Thokas99/veryMAD/tree/v0.5.0) - 2026-08-06

### Added

- Introduced explicit metric settings, annotation and report output, and
  robust scaling for vectors and matrices.
- Added Seurat metadata support, raw-scale thresholds, and unresolved
  flags for zero MAD.

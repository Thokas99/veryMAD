# Changelog

All notable changes to veryMAD are documented here. This changelog follows
[Keep a Changelog 1.1.0](https://keepachangelog.com/en/1.1.0/) and
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Added raw-scale MAD cutoff overlays to `plot_mad_qc_scatter()` when given a
  matching `mad_qc()` report.
- Added a 10,000-cell SPARSim vignette with a composed, batch-aware QC plot.

### Changed

- Reduced point size and opacity to keep dense QC plots readable.

## [0.7.0] - 2026-10-14

### Added

- Added batch-specific and reference-set thresholds to `mad_qc()`.
- Added `plot_mad_qc()` for batch distributions and `plot_mad_qc_scatter()`
  for paired metrics with per-batch Spearman rho and finite-pair counts.

### Changed

- Clarified Seurat metric definitions and QC review guidance.

## [0.6.0] - 2026-10-14

### Added

- Added threshold counts and report-only output to `mad_qc()`.
- Added single-cell examples and a reproducible benchmark.

### Changed

- Set `mad_scale()` to `zero_mad = "na"` by default; `"zero"` uses a unit
  denominator and returns unscaled deviations when centered.
- Reserved `id` and `mad_qc_outlier` as metric names in `mad_qc()`.

## [0.5.1] - 2026-08-06

### Added

- Established `mad_qc()` and `mad_scale()` as the public API.
- Added the combined `mad_qc_outlier` flag and compact report printing.

## [0.5.0] - 2026-08-06

### Added

- Introduced explicit metric settings, annotation and report output, and robust
  scaling for vectors and matrices.
- Added Seurat metadata support, raw-scale thresholds, and unresolved flags for
  zero MAD.

[Unreleased]: https://github.com/Thokas99/veryMAD/compare/v0.7.0...HEAD
[0.7.0]: https://github.com/Thokas99/veryMAD/compare/v0.6.0...v0.7.0
[0.6.0]: https://github.com/Thokas99/veryMAD/compare/v0.5.1...v0.6.0
[0.5.1]: https://github.com/Thokas99/veryMAD/compare/v0.5.0...v0.5.1
[0.5.0]: https://github.com/Thokas99/veryMAD/tree/v0.5.0

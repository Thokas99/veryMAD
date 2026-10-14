# Changelog

## veryMAD (development version)

- Added batch-specific and reference-set thresholds to
  [`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md).
- Added
  [`plot_mad_qc()`](https://thokas99.github.io/veryMAD/reference/plot_mad_qc.md)
  to review metric distributions and flags by batch.
- Clarified Seurat metric definitions and QC review guidance.

## veryMAD 0.6.0

- Set
  [`mad_scale()`](https://thokas99.github.io/veryMAD/reference/mad_scale.md)
  to `zero_mad = "na"` by default; `"zero"` uses a unit denominator and
  returns unscaled deviations when centered.
- Reserved `id` and `mad_qc_outlier` as metric names in
  [`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md).
- Added threshold counts, vectorized combined flags, and report-only
  output.
- Added single-cell examples, regression tests, an optional `scuttle`
  comparison, and a reproducible benchmark.

## veryMAD 0.5.1

- Established
  [`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md)
  and
  [`mad_scale()`](https://thokas99.github.io/veryMAD/reference/mad_scale.md)
  as the public API.
- Added the combined `mad_qc_outlier` flag and compact report printing.
- Expanded input and edge-case coverage and updated package guides.

## veryMAD 0.5.0

- Introduced explicit metric settings, annotation and report output, and
  robust scaling for vectors and matrices.
- Added Seurat metadata support, raw-scale thresholds, and unresolved
  flags for zero MAD.
- Added package guides and streamlined internal calculations.

# veryMAD 0.6.0

- `mad_scale()` defaults to `zero_mad = "na"`; `"zero"` uses a unit
  denominator and returns unscaled deviations when centered.
- `mad_qc()` reserves `id` and `mad_qc_outlier` as metric names.
- Threshold reports include used, flagged, and unresolved counts. Combined
  flags are vectorized, and report-only calls skip annotation.
- Added single-cell examples, regression tests, an optional `scuttle`
  comparison, and a reproducible benchmark script.

# veryMAD 0.5.1

- The public API is `mad_qc()` and `mad_scale()`; QC reports include an overall
  `mad_qc_outlier` flag.
- Added compact report printing and expanded input and edge-case coverage.
- Updated the README, guides, and package checks.

# veryMAD 0.5.0

- Centered the API on `mad_qc()` and `mad_scale()`, with annotation or report
  output, explicit metric settings, and robust scaling for vectors and matrices.
- Added Seurat metadata support and set raw-scale thresholds and unresolved
  zero-MAD flags as QC defaults.
- Added focused guides and streamlined internal calculations.

# Changelog

## veryMAD 0.3.0

### Added

- `summarize_mad_qc()` for deterministic metric-level and observation-level
  summaries of flags already present in a tidy QC report.
- `plot_mad_qc()` for explicit distribution and observation-index views with
  grouped, finite threshold indicators through optional ggplot2.
- Safe `Matrix::dgCMatrix` support in `row_mad()`, `col_mad()`, and
  `top_mad_features()` through optional sparseMatrixStats methods.

### Changed

- Reframed `mad_qc()` documentation around assay-agnostic observation metadata:
  bulk libraries, single cells, or spatial observations remain ungrouped by
  default, while explicit `group_by` changes the statistical reference distribution.
- `mad_qc()` reports now retain calculation-group and compact observation-metadata
  attributes so plots can distinguish saved threshold groups from optional visual faceting.
- `plot_mad_qc()` now uses `facet_by` only for layout and reads calculation groups
  from the report; ungrouped reports never acquire sample or batch thresholds.
- `top_mad_features()` now explicitly preserves original feature order when
  MAD values are tied.
- ggplot2, Matrix, and sparseMatrixStats are optional suggested dependencies;
  dense workflows and package installation do not require them.
- `robust_scale()` continues to reject sparse matrices with a clearer message:
  median centering can destroy sparsity, and veryMAD never converts implicitly
  to a dense matrix.

## veryMAD 0.2.0

### Added

- `mad_limits()` as the shared, auditable threshold foundation.
- `col_mad()` for strictly column-oriented matrix MADs.
- `annotate_mad_qc()` for index-safe metadata annotation.
- `mad_qc_seurat()` for explicit SeuratObject reporting and annotation.
- `top_mad_features()` for deliberately limited raw marginal MAD ranking.
- A focused QC vignette and Linux R-release package-check workflow.

### Changed

- `mad_qc()` is now an ordinary data-frame function with long, stable,
  group-auditable results; multiple grouping columns and missing group labels
  are supported without dropping observations.
- Scores, thresholds, outliers, winsorization, robust scaling, and grouped QC
  now share `zero_mad = "zero"`, `"na"`, or `"error"` behavior.
- Transformations are validated: `log10` requires positive values and `log10p`
  requires nonnegative values.
- `row_mad()` is row-only, `col_mad()` is column-only, and `robust_scale()` uses
  `margin = "rows"` or `"columns"`.
- Runtime dependencies were reduced to `matrixStats`; SeuratObject is optional.

### Breaking changes

- `is_outlier()` was removed; use `is_mad_outlier()` with `nmads` and `direction`.
- `winsorize()` was removed; use `winsorize_mad()`.
- `select_variable_features()` was replaced by `top_mad_features()`, which ranks
  rows only and documents that raw MAD is not mean-variance-aware.
- `mad_qc(Seurat object)` was removed; use `mad_qc_seurat()` explicitly.
- S7 dispatch, dynamic method registration, `.onLoad()` hooks, progress bars,
  `verbose`, numeric matrix margins, and direct Seurat slot manipulation were removed.

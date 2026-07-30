# Changelog

## veryMAD 0.4.0

### Added

- Added `mad_qc_bulk()` and `mad_qc_sc()` as explicit bulk and single-cell
  observation-level wrappers backed by compact `flags`, `thresholds`, and
  `settings` results.
- Added default `log1p` calculation, partial `"none"` overrides, raw-scale
  inverse thresholds, concise `cli` summaries, deterministic playful grades,
  and an optional interactive ASCII startup banner.
- Added coverage for invalid metrics, transformations, matrices, grouped and
  zero-MAD calculations, three-valued flags, collisions, Seurat preservation,
  CLI output, and the banner helper.

### Changed

- Bumped the package version to 0.4.0 and added `cli` to Imports.
- All wrapper metrics are mandatory and explicit; no input class, column name,
  regular expression, or numeric type triggers metric selection.
- Single-cell annotation adds only `mad_qc_outlier`; metric flags remain in the
  compact report. Neither wrapper filters observations.
- `mad_qc_seurat()` now requires explicit metrics and remains as a documented
  compatibility path for the v0.3.2 long-report/per-metric annotation behavior.
- Rewrote the README and vignette around the two-wrapper interface and explicit
  scientific interpretation.

## veryMAD 0.3.2

### Added

- Added `mad_scale()` as the single focused dense-matrix operation, using
  matrixStats row and column medians and MADs for expression heatmap scaling.
- Added realistic deterministic bulk RNA-seq and single-cell metadata examples.

### Changed

- Reports from `mad_limits()`, `mad_qc()`, and `summarize_mad_qc()` are valid
  tibbles; `mad_qc` reports retain calculation metadata for auditing.
- Inactive directional limits are `NA_real_` rather than `Inf` or `-Inf`, and
  outlier evaluation is explicitly direction-aware.
- Observation summaries and annotations use three-valued overall logic:
  `TRUE` for any failure, `FALSE` only when everything passed, and `NA` when
  nothing failed but at least one metric was not evaluated.
- `plot_mad_qc()` uses green, blue, and orange status colors; missing values are
  visible as blue panel-floor triangles; threshold types use dark labeled lines.
- `mad_qc_seurat()` defaults to report mode and both data-frame and Seurat
  annotation protect existing flag columns unless `overwrite = TRUE`.
- Documentation now presents one assay-agnostic metadata-table model with
  ungrouped QC as the universal default and optional deliberate grouping.

### Removed

- Removed sparse-matrix support and the Matrix and sparseMatrixStats suggestions.
- Removed `row_mad()`, `col_mad()`, and raw-MAD feature ranking through
  `top_mad_features()`.
- Replaced the broad `robust_scale()` API with focused `mad_scale()`.
- Removed `winsorize_mad()`; veryMAD flags observations rather than modifying
  measured values.

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
- Transformations are validated: `log10` requires positive values and `log1p`
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

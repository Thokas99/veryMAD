# Changelog

## veryMAD 0.5.1

### API polish

- Reduced the exported API to exactly
  [`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md)
  and
  [`mad_scale()`](https://thokas99.github.io/veryMAD/reference/mad_scale.md).
- Removed the `mad_z_score()` alias.
- Removed the `combine` argument.
- [`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md)
  now always returns the overall `mad_qc_outlier` flag.

### Reliability

- Added compact printing for `verymad_qc` reports.
- Hardened empty-input, missing-value, identifier, overwrite, and
  zero-MAD handling.
- Expanded tests for the complete public API and important edge cases.
- Expanded continuous integration across Linux, macOS, Windows, and
  multiple R versions.

### Documentation

- Updated the README and vignettes for the final two-function interface.
- Simplified the pkgdown reference index.
- Kept `pak::pak("Thokas99/veryMAD")` as the primary installation
  command.

## veryMAD 0.5.0

### Breaking changes

- Reduced the public API to
  [`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md)
  and
  [`mad_scale()`](https://thokas99.github.io/veryMAD/reference/mad_scale.md).
- Removed grouped QC calculations and separate modality workflows.
- Removed redundant per-result counts and qualitative interpretation
  fields.
- Changed the default QC transformation to `"none"`.
- Changed the default zero-MAD policy to `"na"`.

### New interface

- [`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md)
  provides annotation and compact report outputs.
- QC metrics, directions, and transformations remain explicit.
- Seurat support delegates to the same canonical engine.
- [`mad_scale()`](https://thokas99.github.io/veryMAD/reference/mad_scale.md)
  supports robust scaling for vectors and matrices.

### Documentation

- Simplified the README.
- Added focused vignettes for QC, interpretation, and MAD scaling.
- Reorganized the pkgdown website around the two-function API.

### Internal changes

- Consolidated QC calculations into one internal engine.
- Converted low-level statistical functions to internal helpers.
- Removed duplicated wrappers and obsolete implementation paths.
- Simplified dependencies and tests.

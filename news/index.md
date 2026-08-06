# Changelog

## veryMAD 0.5.0

### Breaking changes

- Reduced the public API to
  [`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md)
  and
  [`mad_scale()`](https://thokas99.github.io/veryMAD/reference/mad_scale.md).
- Added
  [`mad_z_score()`](https://thokas99.github.io/veryMAD/reference/mad_scale.md)
  as an alias for robust MAD scaling.
- Removed grouped QC calculations and separate modality workflows.
- Removed redundant per-result counts and qualitative QC grades.
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

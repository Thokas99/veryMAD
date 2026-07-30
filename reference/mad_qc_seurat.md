# Run explicit MAD QC on Seurat cell metadata

`mad_qc_seurat()` preserves the v0.3.2 long report and per-metric
annotation contract. New code should use
[`mad_qc_sc()`](https://thokas99.github.io/veryMAD/reference/mad_qc_sc.md)
for compact reports and one-column annotation. Neither function filters
cells.

## Usage

``` r
mad_qc_seurat(
  object,
  metrics = NULL,
  nmads = 3,
  transform = NULL,
  group_by = NULL,
  zero_mad = c("zero", "na", "error"),
  action = c("report", "annotate"),
  overwrite = FALSE,
  verbose = FALSE
)
```

## Arguments

- object:

  A Seurat object.

- metrics:

  A named character vector mapping metadata columns to `"lower"`,
  `"upper"`, or `"both"`.

- nmads:

  Number of MADs from the median used to define thresholds.

- transform:

  Optional named character vector of explicit per-metric transformations
  passed to
  [`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md).

- group_by:

  Optional metadata column used to estimate thresholds within groups.

- zero_mad:

  Compatibility policy for groups with zero MAD: `"zero"`, `"na"`, or
  `"error"`.

- action:

  Return a QC `"report"` or `"annotate"` the object metadata.

- overwrite:

  Allow annotation columns to replace existing columns.

- verbose:

  Print the concise QC summary?

## Value

For `action = "report"`, a `mad_qc` data frame. For
`action = "annotate"`, the input Seurat object with metadata flag
columns added.

## Details

Use explicit transformations for count-like cell metrics when needed.
For example, `nCount_RNA` and `nFeature_RNA` are commonly right-skewed
and can be thresholded on `log1p`, while bounded percentages such as
`percent.mt` usually remain on the raw scale. Upper-tail count or
feature flags are inspection warnings, not doublet calls.

## Examples

``` r
metrics <- c(
  nCount_RNA = "both",
  nFeature_RNA = "both",
  percent.mt = "upper"
)
transform <- c(
  nCount_RNA = "log1p",
  nFeature_RNA = "log1p"
)
# report <- mad_qc_seurat(seurat_object, metrics, transform = transform)
# seurat_object <- mad_qc_seurat(
#   seurat_object, metrics, transform = transform, action = "annotate"
# )
```

# veryMAD

![veryMAD R package logo](reference/figures/veryMAD-logo.svg)

`veryMAD` 0.4.0 provides small MAD vector primitives plus two explicit
observation-level QC wrappers. It never guesses metrics and never
filters samples, cells, libraries, or spots.

## Contents

1.  [Install](#install)
2.  [Public API](#public-api)
3.  [Bulk QC](#bulk-qc)
4.  [Single-cell QC](#single-cell-qc)
5.  [Transformations and groups](#transformations-and-groups)
6.  [Interpretation](#interpretation)
7.  [Compatibility](#compatibility)

## Install

``` r

pak::pak("Thokas99/veryMAD")
```

Set `options(veryMAD.quiet = TRUE)` to suppress the interactive startup
banner.

## Public API

- [`mad_score()`](https://thokas99.github.io/veryMAD/reference/mad_score.md),
  [`mad_limits()`](https://thokas99.github.io/veryMAD/reference/mad_limits.md),
  [`is_mad_outlier()`](https://thokas99.github.io/veryMAD/reference/is_mad_outlier.md),
  and
  [`mad_scale()`](https://thokas99.github.io/veryMAD/reference/mad_scale.md)
  are focused vector/matrix primitives.
- [`mad_qc_bulk()`](https://thokas99.github.io/veryMAD/reference/mad_qc_bulk.md)
  annotates bulk or other observation metadata with one flag per
  selected metric.
- [`mad_qc_sc()`](https://thokas99.github.io/veryMAD/reference/mad_qc_sc.md)
  reports on cell metadata or adds one combined flag.

Rows are observations. Columns may contain metadata and candidate QC
metrics, but `metrics` must select every metric explicitly:

``` r

metrics <- c(library_size = "lower", mapping_rate = "lower",
             duplication_rate = "upper")
```

Directions are `"lower"`, `"upper"`, or `"both"`. An ordinary
gene-by-sample expression matrix is not automatically a QC-metric table.

## Bulk QC

``` r

bulk <- data.frame(
  library_size = c(1e6, 30e6, 32e6, 31e6, 33e6),
  mapping_rate = c(.50, .92, .94, .93, .95),
  batch = c("a", "a", "a", "b", "b"),
  row.names = paste0("sample", 1:5)
)

bulk_flagged <- mad_qc_bulk(
  bulk,
  metrics = c(library_size = "lower", mapping_rate = "lower"),
  transform = c(mapping_rate = "none")
)

report <- attr(bulk_flagged, "mad_qc")
report$flags
report$thresholds
```

Only `library_size_mad_outlier` and `mapping_rate_mad_outlier` are
appended. The original values, row order, and identifiers are retained.
A numeric matrix works the same way and is returned as an annotated data
frame:

``` r

mad_qc_bulk(as.matrix(bulk[c("library_size", "mapping_rate")]),
            c(library_size = "lower", mapping_rate = "lower"),
            transform = c(mapping_rate = "none"), verbose = FALSE)
```

Inspect flags before making a separate, documented filtering decision:

``` r

bulk_flagged[bulk_flagged$library_size_mad_outlier %in% TRUE, ]
```

## Single-cell QC

``` r

cells <- data.frame(
  nCount_RNA = c(NA, 500, 520, 540, 5000),
  nFeature_RNA = c(NA, 200, 210, 220, 1000),
  percent.mt = c(NA, 3, 4, 5, 30),
  row.names = paste0("cell", 1:5)
)

sc_report <- mad_qc_sc(
  cells,
  metrics = c(nCount_RNA = "lower", nFeature_RNA = "lower", percent.mt = "upper"),
  transform = c(percent.mt = "none"),
  action = "report",
  verbose = FALSE
)

annotated_cells <- mad_qc_sc(
  cells,
  metrics = c(nCount_RNA = "lower", nFeature_RNA = "lower", percent.mt = "upper"),
  transform = c(percent.mt = "none")
)
```

Annotation adds only `mad_qc_outlier`: `TRUE` means at least one
evaluated metric failed, `FALSE` means every metric was evaluated and
passed, and `NA` means none failed but at least one could not be
evaluated. Metric flags remain in `sc_report$flags`; thresholds are in
`sc_report$thresholds`.

Seurat support is optional through `SeuratObject`:

``` r

seurat_object <- mad_qc_sc(
  seurat_object,
  metrics = c(nCount_RNA = "lower", nFeature_RNA = "lower", percent.mt = "upper"),
  transform = c(percent.mt = "none"),
  action = "annotate"
)
```

Cells, assays, reductions, identities, and existing metadata are
preserved.

## Transformations and groups

The wrapper default is `transform = "log1p"` for every selected metric.
A named vector is a partial override, so this leaves `library_size` on
`log1p`:

``` r

transform <- c(mapping_rate = "none", duplication_rate = "none")
```

Use scalar `"none"` to disable all transformations. `"identity"` remains
an alias for `"none"`, and `"log10"` remains supported. Transformations
affect only threshold estimation and calls; raw measurements are
unchanged. Bounded percentages, proportions, and rates may be easier to
interpret with `"none"`, but veryMAD does not decide that for you.

Grouping is also explicit:

``` r

grouped <- mad_qc_bulk(bulk,
  c(library_size = "lower", mapping_rate = "lower"),
  transform = c(mapping_rate = "none"), group_by = "batch", verbose = FALSE)
attr(grouped, "mad_qc")$thresholds
```

The compact result stores flags once, one threshold row per
metric/group, and calculation settings; it does not repeat thresholds
for every observation.

## Interpretation

The CLI grade is playful presentation based on the proportion of unique
observations flagged: 0%, \>0–5%, \>5–15%, and \>15%. It is not a
biological conclusion, acceptance criterion, claim of laboratory
failure, or filtering recommendation.

MAD thresholds are adaptive heuristics, not universal biological
cutoffs. Upper-tail RNA count flags are not automatically doublet calls;
doublet detection is outside veryMAD. Current Bioconductor RNA-QC tools
similarly use robust, metric-specific outlier rules, while veryMAD
deliberately leaves metric, direction, transformation, grouping, and
downstream action under user control. See the [scrapper reference
manual](https://bioconductor.org/packages/release/bioc/manuals/scrapper/man/scrapper.pdf).

## Compatibility

[`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md),
[`annotate_mad_qc()`](https://thokas99.github.io/veryMAD/reference/annotate_mad_qc.md),
[`summarize_mad_qc()`](https://thokas99.github.io/veryMAD/reference/summarize_mad_qc.md),
and
[`plot_mad_qc()`](https://thokas99.github.io/veryMAD/reference/plot_mad_qc.md)
keep the v0.3.2 long-report workflow.
[`mad_qc_seurat()`](https://thokas99.github.io/veryMAD/reference/mad_qc_seurat.md)
remains available, requires explicit metrics, and retains its historical
report/per-metric annotation behavior; new code should use
[`mad_qc_sc()`](https://thokas99.github.io/veryMAD/reference/mad_qc_sc.md)
for the compact report and one-column annotation contract.

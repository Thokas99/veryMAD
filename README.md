<p align="center">
  <img src="man/figures/veryMAD-logo.svg" alt="veryMAD logo" width="360">
</p>

# veryMAD

`veryMAD` is a lightweight, assay-agnostic R toolkit for median absolute
deviation (MAD) quality control. It provides vector methods, dense and sparse
matrix MAD operations, tidy observation-level QC reports, summaries, plots, and
an optional Seurat bridge. It never filters observations automatically.

## Table of contents

1. [What veryMAD does](#what-verymad-does)
2. [Observation-level QC data model](#observation-level-qc-data-model)
3. [Vector-level MAD functions](#vector-level-mad-functions)
4. [Bulk RNA-seq sample QC](#bulk-rna-seq-sample-qc)
5. [Single-cell QC for one object](#single-cell-qc-for-one-object)
6. [Optional grouped QC for merged objects](#optional-grouped-qc-for-merged-objects)
7. [Summaries and plots](#summaries-and-plots)
8. [Seurat integration](#seurat-integration)
9. [Sparse matrix operations](#sparse-matrix-operations)
10. [Zero-MAD behavior](#zero-mad-behavior)
11. [Limitations](#limitations)

## What veryMAD does

Install the development version from GitHub with pak:

```r
install.packages("pak")
pak::pak("Thokas99/veryMAD")
```

The core workflow is:

```text
observation metadata → mad_qc() → inspect → summarize → plot → annotate
```

The result is an ordinary long data frame that can be inspected, saved, joined,
and audited without tidyverse dependencies.

## Observation-level QC data model

`mad_qc()` accepts an observation-level QC metadata table:

```text
rows    = observations
columns = QC metrics and optional grouping variables
```

The observational unit depends on the assay:

```text
bulk RNA-seq:     one row = one sample or sequencing library
single-cell RNA:  one row = one cell
spatial data:     one row = one spot or cell
```

Do not pass a gene-expression count matrix to `mad_qc()`. For bulk RNA-seq,
first calculate sample-level metrics such as library size, detected genes,
mapping rate, duplication rate, mitochondrial read proportion, and ribosomal
read proportion. Pass that sample metadata table to `mad_qc()`.

The canonical interface remains:

```r
mad_qc(data, metrics, nmads = 3, group_by = NULL, ...)
```

`group_by` is a statistical operation, not an organizational convenience:

- `group_by = NULL` calculates one global median, MAD, lower threshold, upper
  threshold, and outlier status per metric.
- `group_by = "sample_id"` calculates those quantities separately within each
  sample.
- `group_by = c("sample_id", "batch")` calculates them within each
  sample–batch combination.

> Small groups produce unstable MAD estimates. Grouped QC should not be used
> merely because a grouping column exists.

## Vector-level MAD functions

```r
library(veryMAD)
x <- c(sample_a = 10, sample_b = 11, sample_c = 12, sample_d = 80)

mad_score(x)
mad_limits(x, nmads = 2, direction = "upper")
is_mad_outlier(x, nmads = 2, direction = "upper")
winsorize_mad(x, nmads = 2, direction = "upper")
```

`direction = "lower"` checks only small values, `"upper"` only large values,
and `"both"` both tails.

## Bulk RNA-seq sample QC

This deterministic documentation dataset contains 250 sequencing libraries,
including low-library, low-complexity, low-mapping, high-duplication,
multi-metric failures, and missing values. Each row is one library.

```r
bulk_metadata <- veryMAD:::.simulate_bulk_qc_metadata(n = 250, seed = 123)

bulk_qc_report <- bulk_metadata |>
  mad_qc(
    metrics = c(
      library_size = "lower",
      detected_genes = "lower",
      mapping_rate = "lower",
      duplication_rate = "upper",
      percent_mitochondrial = "upper"
    ),
    transform = c(
      library_size = "log10p",
      detected_genes = "log10p"
    )
  )

bulk_qc_report |>
  summarize_mad_qc(level = "metric")

bulk_sample_summary <- bulk_qc_report |>
  summarize_mad_qc(level = "observation")

plot_mad_qc(bulk_qc_report, type = "distribution")
plot_mad_qc(bulk_qc_report, type = "index")
```

The 250 libraries form one reference distribution per metric. `condition` is
biological metadata and does not automatically create QC groups or plot
divisions. Grouping bulk samples by condition may hide condition-associated
technical or biological shifts; it should be an explicit analytical choice.

## Single-cell QC for one object

This separate deterministic dataset contains 600 cells from one sample, with
low-count, low-feature, high-mitochondrial, high-doublet, multi-metric, and
missing observations. It deliberately has no fake sample grouping column.

```r
cell_metadata <- veryMAD:::.simulate_single_cell_qc_metadata(n = 600, seed = 123)

cell_qc_report <- cell_metadata |>
  mad_qc(
    metrics = c(
      nCount_RNA = "lower",
      nFeature_RNA = "lower",
      percent.mt = "upper",
      doublet_score = "upper"
    ),
    transform = c(
      nCount_RNA = "log10p",
      nFeature_RNA = "log10p"
    )
  )
```

All 600 cells form one reference distribution per metric. The equivalent
Seurat report is also ungrouped by default:

```r
cell_qc_report <- mad_qc_seurat(
  sample_object,
  metrics = c(
    nCount_RNA = "lower",
    nFeature_RNA = "lower",
    percent.mt = "upper"
  ),
  transform = c(
    nCount_RNA = "log10p",
    nFeature_RNA = "log10p"
  ),
  action = "report"
)
```

## Optional grouped QC for merged objects

Grouped QC can be appropriate for a merged object when samples were processed
independently, depth differs strongly, mitochondrial distributions differ for
technical reasons, every sample has enough cells, and a global threshold would
unfairly remove cells from one sample.

```r
merged_qc_report <- mad_qc_seurat(
  merged_object,
  metrics = c(
    nCount_RNA = "lower",
    nFeature_RNA = "lower",
    percent.mt = "upper"
  ),
  transform = c(
    nCount_RNA = "log10p",
    nFeature_RNA = "log10p"
  ),
  group_by = "orig.ident",
  action = "report"
)
```

Grouping may be inappropriate when groups are small, group differences are
important QC signals, grouping hides a globally poor sample, biological groups
are mistaken for technical batches, or thresholds become incomparable. Merged
Seurat objects do not automatically require sample-specific thresholds.

## Summaries and plots

```r
metric_summary <- bulk_qc_report |>
  summarize_mad_qc(level = "metric")

observation_summary <- bulk_qc_report |>
  summarize_mad_qc(level = "observation")

plot_mad_qc(bulk_qc_report, type = "distribution")
plot_mad_qc(bulk_qc_report, type = "index")
plot_mad_qc(bulk_qc_report, type = "index", facet_by = "condition")
```

Summaries and plots use flags and thresholds already calculated by `mad_qc()`.
They do not recalculate QC. `plot_mad_qc()` reads calculation grouping from the
report attributes only when group-specific thresholds must be shown.

```text
calculation grouping = which observations share a MAD reference distribution
visual faceting      = how observations are arranged in a plot
```

`facet_by` changes only plot layout. It never changes thresholds. An ungrouped
report has one distribution per metric by default—never one division per sample.

## Seurat integration

```r
report <- mad_qc_seurat(sample_object, action = "report")
sample_object <- mad_qc_seurat(sample_object, action = "annotate")
```

The optional bridge uses SeuratObject accessors and the same data-frame
`mad_qc()` implementation. Annotation adds reversible logical metadata flags;
it never filters cells.

## Sparse matrix operations

```r
sparse_counts <- Matrix::rsparsematrix(250, 40, density = 0.05)

row_mad(sparse_counts)
col_mad(sparse_counts)
top_mad_features(sparse_counts, n = 20)
```

`Matrix::dgCMatrix` inputs use optional `sparseMatrixStats` methods without
dense conversion. Sparse raw-MAD ranking is not highly variable gene selection.
`robust_scale()` rejects sparse matrices because median centering can destroy
sparsity.

## Zero-MAD behavior

The default `zero_mad = "zero"` returns neutral zero scores, no outlier flags,
unchanged winsorized values, and zero robust-scaled values for tied margins.
Use `"na"` when the result should be unknown or `"error"` when zero MAD must
stop the analysis. Missing inputs remain missing.

## Limitations

- MAD thresholds are heuristic QC rules, not universal biological truths.
- Small groups produce unstable thresholds; inspect distributions before filtering.
- Grouping can hide global sample failures or meaningful shifts.
- Missing flags mean “Not evaluated,” not passing or failing.
- Raw MAD feature ranking does not model mean–variance dependence.
- Sparse MAD support is limited to `Matrix::dgCMatrix`; sparse robust scaling is unsupported.
- Plot faceting is visual only and does not define calculation groups.
- Users should inspect QC distributions rather than blindly filtering observations.

See `vignette("veryMAD-qc")` for the complete assay-agnostic workflow and
[CHANGELOG.md](CHANGELOG.md) for release changes.

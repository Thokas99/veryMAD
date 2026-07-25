<p align="center">
  <img src="man/figures/veryMAD-logo.svg" alt="veryMAD logo" width="220">
</p>

# veryMAD

`veryMAD` is a lightweight R toolkit for robust MAD scores, dense
expression-matrix scaling and auditable quality-control flagging across bulk,
single-cell and other observation-level omics metadata.

Install the development version with:

```r
install.packages("pak")
pak::pak("Thokas99/veryMAD")
```

## Table of contents

1. [Package purpose](#package-purpose)
2. [Observation-level QC data model](#observation-level-qc-data-model)
3. [Robust vector z-scores](#robust-vector-z-scores)
4. [Robust matrix scaling for heatmaps](#robust-matrix-scaling-for-heatmaps)
5. [Bulk RNA-seq sample QC](#bulk-rna-seq-sample-qc)
6. [Single-cell QC for one object](#single-cell-qc-for-one-object)
7. [Optional grouped QC for merged objects](#optional-grouped-qc-for-merged-objects)
8. [Tibble reports](#tibble-reports)
9. [Summaries and plots](#summaries-and-plots)
10. [Safe annotation](#safe-annotation)
11. [Seurat reporting and annotation](#seurat-reporting-and-annotation)
12. [Zero-MAD behavior](#zero-mad-behavior)
13. [Limitations](#limitations)
14. [Public API](#public-api)

## Package purpose

The package does three things deliberately: robust vector scoring, dense
matrix scaling for expression heatmaps, and transparent QC flagging from
observation metadata. It never filters observations, changes measurements, or
guesses which assay, metric, sample, condition, or batch a user intended.

The core workflow is:

```text
metadata → mad_qc() → inspect → summarize → plot → annotate
```

## Observation-level QC data model

`mad_qc()` accepts a metadata table in which rows are observations and columns
are numeric QC metrics plus optional grouping metadata.

| Context | One row |
| --- | --- |
| Bulk RNA-seq | One sample or sequencing library |
| Single-cell RNA-seq | One cell |
| Spatial data | One spot or cell |
| Other omics | One sample or observation |

A gene-expression matrix is not observation-level QC metadata. For bulk RNA-seq,
first calculate sample metrics such as library size, detected genes, mapping
rate, duplication rate, and mitochondrial proportion, then pass that table to
`mad_qc()`.

## Robust vector z-scores

`mad_score()` uses the median and scaled median absolute deviation instead of
the mean and standard deviation:

```r
library(veryMAD)

x <- c(10, 11, 11, 12, 12, 13, 100)

comparison <- tibble::tibble(
  value = x,
  classical_z = as.numeric(scale(x)),
  robust_z = mad_score(x)
)

comparison
mad_limits(x, nmads = 3, direction = "upper")
is_mad_outlier(x, nmads = 3, direction = "upper")
```

Classical scaling uses the mean and standard deviation; MAD scaling uses the
median and median absolute deviation. Extreme values therefore influence the
robust center and scale less strongly. A large robust score indicates distance
from the robust center, not automatic proof of technical failure.

## Robust matrix scaling for heatmaps

`mad_scale()` is the package's only public matrix operation. It handles numeric
vectors and ordinary dense numeric matrices. For an expression matrix with
genes in rows and samples in columns:

```r
normalized_expression <- matrix(
  log2(c(12, 18, 14, 80, 75, 90, 30, 25, 28, 8, 9, 7) + 1),
  nrow = 3,
  byrow = TRUE,
  dimnames = list(c("GeneA", "GeneB", "GeneC"), paste0("Sample_", 1:4))
)
selected_genes <- rownames(normalized_expression)

heatmap_matrix <- mad_scale(
  normalized_expression[selected_genes, , drop = FALSE],
  margin = "rows"
)

if (requireNamespace("pheatmap", quietly = TRUE)) {
  pheatmap::pheatmap(heatmap_matrix, scale = "none")
}
```

Use `scale = "none"` because `mad_scale()` has already standardized each row.
For bulk RNA-seq, use logCPM, VST, rlog, or another appropriate normalized and
transformed matrix—not unnormalized raw counts. For example, given an existing
edgeR `DGEList` called `dge`:

```r
# logcpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)
# heatmap_matrix <- mad_scale(logcpm[selected_genes, , drop = FALSE], margin = "rows")
```

Very small nonzero MAD values can produce extreme scores. `mad_scale()` does
not clip them. For visualization only, an explicit display range can be made:

```r
display_matrix <- pmax(pmin(heatmap_matrix, 4), -4)
```

## Bulk RNA-seq sample QC

This deterministic simulation represents 150 sequencing libraries. Its ranges
are illustrative, not acceptance thresholds.

```r
bulk_metadata <- veryMAD:::.simulate_bulk_qc_metadata(n = 150, seed = 123)

bulk_report <- bulk_metadata |>
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

head(bulk_report)
```

The full collection of libraries forms one reference distribution per metric.
Condition and batch are not automatic QC groups. Grouping bulk samples by a
biological condition can hide condition-associated technical or biological
shifts and should be an explicit analytical choice.

## Single-cell QC for one object

This separate simulation represents 1,000 cells from one object and contains no
redundant sample grouping column:

```r
cell_metadata <- veryMAD:::.simulate_single_cell_qc_metadata(n = 1000, seed = 123)

cell_report <- cell_metadata |>
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

head(cell_report)
```

All cells in this individual object form the reference distribution.

## Optional grouped QC for merged objects

Grouping changes the median, MAD, thresholds, and outlier calls. It is not a
plotting preference. For a merged object, sample-specific thresholds may be
reasonable when samples were processed independently, technical distributions
differ strongly, and every sample contains enough cells:

```r
# Given an existing merged Seurat object called `merged_object`:
# merged_report <- mad_qc_seurat(
#   merged_object,
#   metrics = c(
#     nCount_RNA = "lower",
#     nFeature_RNA = "lower",
#     percent.mt = "upper"
#   ),
#   transform = c(nCount_RNA = "log10p", nFeature_RNA = "log10p"),
#   group_by = "orig.ident",
#   action = "report"
# )
```

Do not group merely because a column exists. Grouping can hide a globally poor
sample, make thresholds incomparable, or mistake biology for a technical batch.
Small groups produce unstable medians and MAD estimates.

## Tibble reports

`mad_limits()`, `mad_qc()`, and `summarize_mad_qc()` return compact tibbles.
A QC report contains one observation–metric pair per row, including raw and
transformed values, thresholds, direction, and the logical flag. Missing flags
mean “Not evaluated.” Reports retain calculation metadata in attributes for
auditing and plotting.

## Summaries and plots

```r
bulk_report |>
  summarize_mad_qc(level = "metric")

bulk_sample_summary <- bulk_report |>
  summarize_mad_qc(level = "observation")

head(bulk_sample_summary)

if (requireNamespace("ggplot2", quietly = TRUE)) {
  plot_mad_qc(bulk_report, type = "distribution")
  plot_mad_qc(bulk_report, type = "index")
}
```

Pass, Not evaluated, and Outlier are green, blue, and orange. A blue triangle at
the panel floor marks a missing value with no natural y-coordinate. Dark dashed
and dot-dashed lines show lower and upper thresholds. `facet_by` changes only
the plot layout; it never changes the calculation groups.

## Safe annotation

```r
annotated_bulk <- annotate_mad_qc(bulk_metadata, bulk_report)
head(annotated_bulk)
```

Annotation preserves the input rows and measurements and adds one logical flag
per metric plus an overall flag. Overall status is `TRUE` for any failure,
`FALSE` only when every metric was evaluated and passed, and `NA` when nothing
failed but at least one metric was not evaluated. Existing flag columns are
protected unless `overwrite = TRUE` is explicit.

## Seurat reporting and annotation

Only the optional `SeuratObject` package is required. Report mode is the safe
default:

```r
# Given an existing Seurat object called `sample_object`:
# cell_report <- mad_qc_seurat(
#   sample_object,
#   metrics = c(
#     nCount_RNA = "lower",
#     nFeature_RNA = "lower",
#     percent.mt = "upper"
#   ),
#   transform = c(nCount_RNA = "log10p", nFeature_RNA = "log10p"),
#   action = "report"
# )
```

Annotation returns a modified object; it does not mutate the variable in place:

```r
# annotated_object <- mad_qc_seurat(sample_object, action = "annotate")
# sample_object <- mad_qc_seurat(sample_object, action = "annotate") # reassign explicitly
```

Only intended metadata flags are added through `SeuratObject::AddMetaData()`.
Cells are never filtered.

## Zero-MAD behavior

For vectors and matrix margins with zero MAD:

- `zero_mad = "zero"` returns neutral zero scores or `FALSE` flags.
- `zero_mad = "na"` returns missing results where scaling or evaluation is undefined.
- `zero_mad = "error"` stops with a clear message.

Missing input positions remain missing. A margin containing only missing values
is not evaluable. Very small nonzero MAD values are not treated as zero.

## Limitations

- MAD thresholds are heuristic QC rules, not universal biological truths.
- Small groups produce unstable medians, MADs, and thresholds.
- Simulation ranges are illustrative rather than recommended cutoffs.
- `mad_scale()` supports ordinary dense matrices only and can produce extreme
  scores when a nonzero MAD is very small.
- No automatic assay detection, grouping, clipping, filtering, or winsorization
  is performed.
- Users should inspect QC distributions and context before excluding observations.

## Public API

| Function | Purpose | Output |
| --- | --- | --- |
| `mad_score()` | Robust vector z-scores | Numeric vector |
| `mad_scale()` | Robust vector or dense-matrix scaling | Numeric vector or matrix |
| `mad_limits()` | Directional MAD thresholds | One-row tibble |
| `is_mad_outlier()` | Direction-aware outlier flags | Logical vector |
| `mad_qc()` | Observation-level metadata QC | `mad_qc` tibble |
| `summarize_mad_qc()` | Metric or observation summaries | Tibble |
| `plot_mad_qc()` | Distribution or index plots | ggplot object |
| `annotate_mad_qc()` | Append aligned QC flags | Original table class |
| `mad_qc_seurat()` | Report or annotate Seurat metadata | QC tibble or Seurat object |

`winsorize_mad()` was removed. `veryMAD` focuses on scoring, reporting, and
flagging observations rather than modifying measured values.

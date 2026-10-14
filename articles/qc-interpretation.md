# Interpreting MAD QC flags

MAD flags identify selected statistical tails for review. Choose a
direction for each metric:

- `lower` flags unusually small values;
- `upper` flags unusually large values;
- `both` flags either tail.

## Read a flag in context

The example uses common sequencing metrics. `pct_mito` is in percentage
points: 8 denotes 8%.

``` r

qc_metadata <- data.frame(
  sample = paste0("sample_", 1:7),
  library_size = c(2.4e6, 2.5e6, 2.6e6, 2.7e6, 2.8e6, 0.5e6, NA_real_),
  detected_genes = c(12000, 12500, 13100, 12800, 13400, 1600, 12200),
  pct_mito = c(4, 5, 6, 5, 7, 8, 30)
)

qc_result <- mad_qc(
  qc_metadata,
  metrics = c(
    library_size = "lower",
    detected_genes = "lower",
    pct_mito = "upper"
  ),
  output = "report",
  verbose = FALSE
)

qc_result$flags[, c("id", "library_size", "detected_genes",
                    "pct_mito", "mad_qc_outlier")]
#>   id library_size detected_genes pct_mito mad_qc_outlier
#> 1  1        FALSE          FALSE    FALSE          FALSE
#> 2  2        FALSE          FALSE    FALSE          FALSE
#> 3  3        FALSE          FALSE    FALSE          FALSE
#> 4  4        FALSE          FALSE    FALSE          FALSE
#> 5  5        FALSE          FALSE    FALSE          FALSE
#> 6  6         TRUE           TRUE    FALSE           TRUE
#> 7  7           NA          FALSE     TRUE           TRUE
```

The missing library size remains `NA` in its flag. The combined flag is
`NA` when no metric is flagged and at least one metric is unresolved.
Review individual flags with the source metadata.

## Threshold status and transformations

The threshold table reports calculation-scale and raw-scale limits.

``` r

transformed_result <- mad_qc(
  qc_metadata,
  metrics = c(library_size = "lower", pct_mito = "upper"),
  transform = c(library_size = "log1p"),
  output = "report",
  verbose = FALSE
)

transformed_result$thresholds[, c(
  "metric", "direction", "transform", "lower_raw", "upper_raw", "status"
)]
#>         metric direction transform lower_raw upper_raw status
#> 1 library_size     lower     log1p   1961998        NA     ok
#> 2     pct_mito     upper      none        NA   10.4478     ok
```

`min_n` sets the minimum number of observations used for estimation. The
statuses `insufficient_n` and `all_missing` identify metrics without an
estimated threshold.

## Zero MAD and downstream decisions

A zero MAD means the data provide no non-zero robust spread estimate.
Choose the policy for the downstream analysis.

``` r

constant_metric <- data.frame(
  sample = paste0("sample_", 1:6),
  detected_genes = rep(12000, 6)
)

zero_mad_result <- mad_qc(
  constant_metric,
  metrics = c(detected_genes = "lower"),
  zero_mad = "zero",
  output = "report",
  verbose = FALSE
)

zero_mad_result$thresholds[, c("metric", "mad", "lower_raw", "status")]
#>           metric mad lower_raw   status
#> 1 detected_genes   0     12000 zero_mad
```

Use `zero_mad = "na"` to leave flags unresolved or `"error"` to stop.
For stratified analyses, split the data and call
[`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md) on
each subset.

## Single-cell RNA-seq

Common RNA metrics are `nCount_RNA`, `nFeature_RNA`, and `percent.mt`.
This example checks lower count and feature tails, an upper
mitochondrial tail, and applies `log1p` to the count-like metrics.
Seurat’s `percent.mt` is usually on a 0–100 scale: 8 means 8%.

``` r

rna_metadata <- data.frame(
  sample = rep(c("sample_A", "sample_B"), each = 6),
  nCount_RNA = c(820, 910, 970, 1050, 1120, 100, 125, 150, 175, 190, 210, 230),
  nFeature_RNA = c(430, 480, 520, 570, 610, 55, 80, 95, 100, 115, 125, 140),
  percent.mt = c(5, 6, 7, 8, 9, 35, 4, 5, 6, 7, 8, 9)
)

rna_metrics <- c(
  nCount_RNA = "lower",
  nFeature_RNA = "lower",
  percent.mt = "upper"
)
rna_transform <- c(nCount_RNA = "log1p", nFeature_RNA = "log1p")

rna_report <- mad_qc(
  rna_metadata,
  metrics = rna_metrics,
  transform = rna_transform,
  output = "report",
  verbose = FALSE
)
rna_report$thresholds[, c("metric", "direction", "transform", "status")]
#>         metric direction transform status
#> 1   nCount_RNA     lower     log1p     ok
#> 2 nFeature_RNA     lower     log1p     ok
#> 3   percent.mt     upper      none     ok
```

## Single-cell ATAC-seq

ATAC metadata may include `nCount_peaks`, `TSS.enrichment`, and
`nucleosome_signal`. The example selects lower peak-count and TSS tails
and an upper nucleosome-signal tail. Names and useful directions depend
on the pipeline and assay; these are not universal thresholds.

``` r

atac_metadata <- data.frame(
  nCount_peaks = c(1200, 1800, 3200, 6000, 9000, 20000, 2500, 4100),
  TSS.enrichment = c(0.8, 1.3, 3.0, 5.4, 8.0, 15.0, 4.1, 6.2),
  nucleosome_signal = c(3.5, 2.4, 1.8, 1.4, 1.1, 0.9, 1.6, 1.3)
)

atac_report <- mad_qc(
  atac_metadata,
  metrics = c(
    nCount_peaks = "lower",
    TSS.enrichment = "lower",
    nucleosome_signal = "upper"
  ),
  transform = c(nCount_peaks = "log1p"),
  output = "report",
  verbose = FALSE
)
atac_report$thresholds[, c("metric", "direction", "transform", "status")]
#>              metric direction transform status
#> 1      nCount_peaks     lower     log1p     ok
#> 2    TSS.enrichment     lower      none     ok
#> 3 nucleosome_signal     upper      none     ok
```

## Per-sample QC and absolute criteria

Estimate thresholds per sample by splitting metadata and calling
[`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md) on
each subset. Sample-specific thresholds can hide a uniformly low-quality
sample; also review between-sample distributions.

``` r

sample_reports <- lapply(
  split(rna_metadata, rna_metadata$sample),
  function(x) mad_qc(
    x,
    metrics = rna_metrics,
    transform = rna_transform,
    output = "report",
    verbose = FALSE
  )
)
sample_reports$sample_A$thresholds[, c("metric", "lower_raw", "status")]
#>         metric lower_raw status
#> 1   nCount_RNA  542.0431     ok
#> 2 nFeature_RNA  266.8073     ok
#> 3   percent.mt        NA     ok
```

MAD flags identify outliers in the supplied reference distribution.
Combine them with separately chosen absolute criteria when appropriate.
These example cutoffs are illustrative, not general recommendations:

``` r

rna_annotated <- mad_qc(
  rna_metadata,
  metrics = rna_metrics,
  transform = rna_transform,
  verbose = FALSE
)

illustrative_keep <- !rna_annotated$mad_qc_outlier &
  rna_annotated$nCount_RNA >= 500 &
  rna_annotated$percent.mt <= 20
table(illustrative_keep, useNA = "ifany")
#> illustrative_keep
#> FALSE  TRUE 
#>     7     5
```

## Seurat metadata

[`mad_qc()`](https://thokas99.github.io/veryMAD/reference/mad_qc.md)
adds flags to Seurat cell metadata. The example checks cell order
against the original identities.

``` r

if (requireNamespace("SeuratObject", quietly = TRUE)) {
  counts <- matrix(
    seq_len(40) %% 6,
    nrow = 5,
    dimnames = list(paste0("gene_", 1:5), paste0("cell_", 1:8))
  )
  seurat_object <- SeuratObject::CreateSeuratObject(counts = counts)
  seurat_object[["percent.mt"]] <- c(4, 5, 6, 7, 8, 9, 10, 30)
  cell_ids <- colnames(seurat_object)

  seurat_metadata <- seurat_object[[]]
  expected_flags <- mad_qc(
    seurat_metadata,
    metrics = rna_metrics,
    transform = rna_transform,
    verbose = FALSE
  )
  seurat_annotated <- mad_qc(
    seurat_object,
    metrics = rna_metrics,
    transform = rna_transform,
    verbose = FALSE
  )

  stopifnot(identical(rownames(seurat_annotated[[]]), cell_ids))
  stopifnot(identical(
    seurat_annotated[[]]$mad_qc_outlier,
    expected_flags$mad_qc_outlier
  ))
  head(seurat_annotated[[]][, c("nCount_RNA", "nFeature_RNA",
                                "percent.mt", "mad_qc_outlier")])
}
#> Warning: Feature names cannot have underscores ('_'), replacing with dashes
#> ('-')
#> Warning: Data is of class matrix. Coercing to dgCMatrix.
#>        nCount_RNA nFeature_RNA percent.mt mad_qc_outlier
#> cell_1         15            5          4             NA
#> cell_2         10            4          5             NA
#> cell_3         11            4          6             NA
#> cell_4         12            4          7             NA
#> cell_5         13            4          8             NA
#> cell_6         14            4          9             NA
```

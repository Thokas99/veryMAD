# veryMAD

`veryMAD` 0.2.0 is a small, pipe-friendly R toolkit for median absolute
deviation (MAD) scores, thresholds, outliers, robust scaling, and auditable
omics quality control. Its core result is an ordinary data frame: easy to
inspect, save, join, and review.

## Installation

```r
install.packages("pak")
pak::pak("Thokas99/veryMAD")
```

## Vector methods

```r
library(veryMAD)
x <- c(sample_a = 10, sample_b = 11, sample_c = 12, sample_d = 80)

mad_score(x)
mad_limits(x, nmads = 2, direction = "upper")
is_mad_outlier(x, nmads = 2, direction = "upper")
winsorize_mad(x, nmads = 2, direction = "upper")
```

`direction = "lower"` checks only small values, `"upper"` only large values,
and `"both"` both tails. Directional winsorization caps only the selected tail.

## Tidy metadata QC

```r
metadata <- data.frame(
  sample_id = rep(c("S1", "S2"), each = 4),
  nCount_RNA = c(800, 900, 1000, 9000, 700, 850, 950, 8000),
  nFeature_RNA = c(400, 450, 500, 520, 350, 420, 480, 490),
  percent.mt = c(2, 3, 4, 30, 2, 3, 5, 25)
)

qc <- metadata |>
  mad_qc(
    metrics = c(
      nCount_RNA = "lower",
      nFeature_RNA = "lower",
      percent.mt = "upper"
    ),
    nmads = 3,
    group_by = "sample_id",
    transform = c(nCount_RNA = "log10p", nFeature_RNA = "log10p")
  )

qc
annotated <- annotate_mad_qc(metadata, qc)
```

The report keeps the raw and transformed values, thresholds, direction,
grouping columns, stable internal observation index, and flag for every
observation-metric pair. `group_by` accepts one or several string column names;
missing group labels form a group and never drop rows. Available transforms are
`identity`, positive-only `log10`, and nonnegative `log10p` (`log10(x + 1)`).

## SeuratObject integration

```r
report <- mad_qc_seurat(seurat_object, action = "report")
seurat_object <- mad_qc_seurat(seurat_object, action = "annotate")
```

The explicit integration requires only the optional `SeuratObject` package.
Report mode returns the same tidy table; annotation mode adds reversible logical
metadata flags. It never filters cells automatically.

## Zero MAD

All score and threshold functions use `zero_mad = "zero"` by default. A tied
vector or group receives neutral zero scores, no outlier flags, unchanged
winsorized values, and zero robust-scaled values. Use `"na"` when the result
should be unknown or `"error"` when a tied margin must stop the analysis.
Missing input positions remain missing.

## Public API

| Function | Purpose |
| --- | --- |
| `mad_score()` | Robust median/MAD scores. |
| `mad_limits()` | Auditable lower and upper limits. |
| `is_mad_outlier()` | Directional logical flags. |
| `winsorize_mad()` | Directional capping. |
| `row_mad()`, `col_mad()` | Matrix-margin MADs. |
| `robust_scale()` | Row- or column-wise robust scaling. |
| `top_mad_features()` | Rank matrix rows by raw MAD. |
| `mad_qc()` | Long data-frame QC report. |
| `annotate_mad_qc()` | Align report flags back to metadata. |
| `mad_qc_seurat()` | Explicit Seurat report or annotation. |

## Limitations

- MAD thresholds are heuristic QC rules, not universal biological truths.
- Small groups produce unstable thresholds; inspect distributions before filtering.
- Zero-MAD groups require an explicit policy choice.
- Raw MAD feature ranking does not model mean-variance dependence and is not a
  replacement for Seurat or Scanpy highly variable feature selection.
- Sparse matrices are rejected to avoid accidental dense conversion.
- Users should inspect QC distributions rather than blindly filtering cells.

See `vignette("veryMAD-qc")` for the focused report-to-annotation workflow and
[CHANGELOG.md](CHANGELOG.md) for breaking changes.

# Flag observations with explicit MAD thresholds

`mad_qc()` calculates median absolute deviation (MAD) thresholds for
selected numeric metadata columns and returns a long, auditable report.
Rows in `data` are observations: bulk RNA-seq libraries, single cells,
spatial spots, or any other assay-level unit with observation metadata.

## Usage

``` r
mad_qc(
  data,
  metrics,
  nmads = 3,
  group_by = NULL,
  transform = NULL,
  constant = 1.4826,
  na_rm = TRUE,
  zero_mad = c("zero", "na", "error")
)
```

## Arguments

- data:

  A data frame containing one observation per row.

- metrics:

  A named character vector mapping metric names in `data` to one of
  `"lower"`, `"upper"`, or `"both"`.

- nmads:

  Number of MADs from the median used to define thresholds.

- group_by:

  Optional column name used to estimate thresholds within groups.

- transform:

  Optional named character vector mapping metrics to explicit
  transformations. Supported values are `"identity"`, `"log10"`, and
  `"log10p"`. Metrics absent from this vector use identity
  transformation.

- constant:

  Consistency constant passed to MAD calculation.

- na_rm:

  Logical; remove missing values before threshold estimation.

- zero_mad:

  How to handle groups with zero MAD: `"zero"`, `"na"`, or `"error"`.

## Value

A `mad_qc` data frame with one row per observation-metric pair. It
includes `raw_value`, transformed `value`, threshold columns,
`direction`, and `is_outlier`.

## Details

Count-like QC metrics such as library size, total RNA counts, and
detected features are often strongly right-skewed. In those cases,
estimating MAD thresholds after an explicit log transformation such as
`log10p` can give a more useful reference distribution. Bounded
percentages, proportions, and rates are usually kept on their original
scale because a log transform changes the interpretation of already
bounded measurements.

Transformations are explicit and per metric. `transform = NULL` is
equivalent to identity transformation for every metric, and `mad_qc()`
never infers a transformation from a metric name. Thresholds, medians,
MADs, and the `value` column are expressed on the calculation scale. The
`raw_value` column always preserves the exact input measurement scale.

`direction = "both"` flags both tails of the selected metric, but it
does not assign a biological cause. For example, an upper-tail count or
feature flag means unusually high relative to the reference
distribution; it is not a doublet call or an automatic filtering
decision.

## Examples

``` r
bulk_metadata <- data.frame(
  sample_id = paste0("sample_", 1:8),
  library_size = c(2e6, 3e7, 3.2e7, 3.4e7, 3.1e7, 3.3e7, 3.5e7, 3.6e7),
  mapping_rate = c(0.55, 0.92, 0.94, 0.93, 0.95, 0.94, 0.93, 0.92)
)
mad_qc(
  bulk_metadata,
  metrics = c(library_size = "lower", mapping_rate = "lower"),
  transform = c(library_size = "log10p")
)
#> # A tibble: 16 × 11
#>     .obs id    metric        raw_value value median    mad lower upper direction
#>    <int> <chr> <chr>             <dbl> <dbl>  <dbl>  <dbl> <dbl> <dbl> <chr>    
#>  1     1 1     library_size    2   e+6  6.30   7.51 0.0391 7.39     NA lower    
#>  2     2 2     library_size    3   e+7  7.48   7.51 0.0391 7.39     NA lower    
#>  3     3 3     library_size    3.20e+7  7.51   7.51 0.0391 7.39     NA lower    
#>  4     4 4     library_size    3.4 e+7  7.53   7.51 0.0391 7.39     NA lower    
#>  5     5 5     library_size    3.10e+7  7.49   7.51 0.0391 7.39     NA lower    
#>  6     6 6     library_size    3.30e+7  7.52   7.51 0.0391 7.39     NA lower    
#>  7     7 7     library_size    3.5 e+7  7.54   7.51 0.0391 7.39     NA lower    
#>  8     8 8     library_size    3.60e+7  7.56   7.51 0.0391 7.39     NA lower    
#>  9     1 1     mapping_rate    5.5 e-1  0.55   0.93 0.0148 0.886    NA lower    
#> 10     2 2     mapping_rate    9.2 e-1  0.92   0.93 0.0148 0.886    NA lower    
#> 11     3 3     mapping_rate    9.4 e-1  0.94   0.93 0.0148 0.886    NA lower    
#> 12     4 4     mapping_rate    9.3 e-1  0.93   0.93 0.0148 0.886    NA lower    
#> 13     5 5     mapping_rate    9.5 e-1  0.95   0.93 0.0148 0.886    NA lower    
#> 14     6 6     mapping_rate    9.4 e-1  0.94   0.93 0.0148 0.886    NA lower    
#> 15     7 7     mapping_rate    9.3 e-1  0.93   0.93 0.0148 0.886    NA lower    
#> 16     8 8     mapping_rate    9.2 e-1  0.92   0.93 0.0148 0.886    NA lower    
#> # ℹ 1 more variable: is_outlier <lgl>
cell_metadata <- data.frame(
  cell_id = paste0("cell_", 1:8),
  nCount_RNA = c(200, 8000, 8500, 9000, 8700, 9200, 9500, 70000),
  nFeature_RNA = c(150, 2800, 3000, 3100, 2950, 3200, 3300, 8000),
  percent.mt = c(3, 4, 5, 4, 6, 5, 4, 18)
)
mad_qc(
  cell_metadata,
  metrics = c(nCount_RNA = "both", nFeature_RNA = "both", percent.mt = "upper"),
  transform = c(nCount_RNA = "log10p", nFeature_RNA = "log10p")
)
#> # A tibble: 24 × 11
#>     .obs id    metric       raw_value value median    mad lower upper direction
#>    <int> <chr> <chr>            <dbl> <dbl>  <dbl>  <dbl> <dbl> <dbl> <chr>    
#>  1     1 1     nCount_RNA         200  2.30   3.95 0.0358  3.84  4.05 both     
#>  2     2 2     nCount_RNA        8000  3.90   3.95 0.0358  3.84  4.05 both     
#>  3     3 3     nCount_RNA        8500  3.93   3.95 0.0358  3.84  4.05 both     
#>  4     4 4     nCount_RNA        9000  3.95   3.95 0.0358  3.84  4.05 both     
#>  5     5 5     nCount_RNA        8700  3.94   3.95 0.0358  3.84  4.05 both     
#>  6     6 6     nCount_RNA        9200  3.96   3.95 0.0358  3.84  4.05 both     
#>  7     7 7     nCount_RNA        9500  3.98   3.95 0.0358  3.84  4.05 both     
#>  8     8 8     nCount_RNA       70000  4.85   3.95 0.0358  3.84  4.05 both     
#>  9     1 1     nFeature_RNA       150  2.18   3.48 0.0409  3.36  3.61 both     
#> 10     2 2     nFeature_RNA      2800  3.45   3.48 0.0409  3.36  3.61 both     
#> # ℹ 14 more rows
#> # ℹ 1 more variable: is_outlier <lgl>
```

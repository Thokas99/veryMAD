# veryMAD quality-control examples

`veryMAD` calculates explicit MAD-based QC flags on observation-level
metadata. The input measurements are preserved as `raw_value`; optional
transformations only define the calculation scale stored in `value`.

## Choosing a calculation scale

Strongly right-skewed counts, such as library size and detected-feature
counts, are often easier to threshold on a log scale. Bounded
percentages, proportions, and rates usually stay on their original
scale. `veryMAD` never infers this from column names: the transformation
must be visible in the function call.

## Bulk RNA-seq sample QC

``` r

bulk_metadata <- veryMAD:::.simulate_bulk_qc_metadata(n = 240, seed = 123)

bulk_metrics <- c(
  library_size = "lower",
  detected_genes = "lower",
  mapping_rate = "lower",
  assigned_rate = "lower",
  rrna_rate = "upper"
)

bulk_transform <- c(
  library_size = "log10p",
  detected_genes = "log10p"
)

bulk_qc <- mad_qc(
  bulk_metadata,
  metrics = bulk_metrics,
  transform = bulk_transform,
  nmads = 3
)

head(bulk_qc)
#> # A tibble: 6 × 11
#>    .obs id    metric       raw_value value median   mad lower upper direction
#>   <int> <chr> <chr>            <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <chr>    
#> 1     1 1     library_size   3000000  6.48   7.63 0.128  7.24    NA lower    
#> 2     2 2     library_size   4285714  6.63   7.63 0.128  7.24    NA lower    
#> 3     3 3     library_size   5571429  6.75   7.63 0.128  7.24    NA lower    
#> 4     4 4     library_size   6857143  6.84   7.63 0.128  7.24    NA lower    
#> 5     5 5     library_size   8142857  6.91   7.63 0.128  7.24    NA lower    
#> 6     6 6     library_size   9428571  6.97   7.63 0.128  7.24    NA lower    
#> # ℹ 1 more variable: is_outlier <lgl>
summarize_mad_qc(bulk_qc, level = "metric")
#> # A tibble: 5 × 7
#>   metric      n_observations n_evaluated n_missing n_outliers outlier_proportion
#>   <chr>                <int>       <int>     <int>      <int>              <dbl>
#> 1 library_si…            240         240         0          8             0.0333
#> 2 detected_g…            240         240         0         11             0.0458
#> 3 mapping_ra…            240         240         0          8             0.0333
#> 4 assigned_r…            240         240         0          8             0.0333
#> 5 rrna_rate              240         240         0         18             0.075 
#> # ℹ 1 more variable: direction <chr>
```

The simulated ranges are illustrative and are not acceptance cutoffs.
Duplication rate and mitochondrial proportion are not used as universal
bulk defaults here: both can be informative in specific workflows, but
their interpretation depends on protocol and biology.

## Single-cell metadata QC

``` r

cell_metadata <- veryMAD:::.simulate_single_cell_qc_metadata(n = 1200, seed = 123)

sc_metrics <- c(
  nCount_RNA = "both",
  nFeature_RNA = "both",
  percent.mt = "upper"
)

sc_transform <- c(
  nCount_RNA = "log10p",
  nFeature_RNA = "log10p"
)

sc_qc <- mad_qc(
  cell_metadata,
  metrics = sc_metrics,
  transform = sc_transform,
  nmads = 3
)

head(sc_qc)
#> # A tibble: 6 × 11
#>    .obs id    metric     raw_value value median   mad lower upper direction
#>   <int> <chr> <chr>          <dbl> <dbl>  <dbl> <dbl> <dbl> <dbl> <chr>    
#> 1     1 1     nCount_RNA       150  2.18   3.95 0.237  3.24  4.67 both     
#> 2     2 2     nCount_RNA       189  2.28   3.95 0.237  3.24  4.67 both     
#> 3     3 3     nCount_RNA       229  2.36   3.95 0.237  3.24  4.67 both     
#> 4     4 4     nCount_RNA       268  2.43   3.95 0.237  3.24  4.67 both     
#> 5     5 5     nCount_RNA       308  2.49   3.95 0.237  3.24  4.67 both     
#> 6     6 6     nCount_RNA       347  2.54   3.95 0.237  3.24  4.67 both     
#> # ℹ 1 more variable: is_outlier <lgl>
summarize_mad_qc(sc_qc, level = "metric")
#> # A tibble: 3 × 7
#>   metric      n_observations n_evaluated n_missing n_outliers outlier_proportion
#>   <chr>                <int>       <int>     <int>      <int>              <dbl>
#> 1 nCount_RNA            1200        1200         0         37             0.0308
#> 2 nFeature_R…           1200        1200         0         43             0.0358
#> 3 percent.mt            1200        1200         0         61             0.0508
#> # ℹ 1 more variable: direction <chr>
```

The Bioconductor single-cell workflow commonly uses lower-tail filtering
for library size and detected features. The example uses two-sided count
and feature flags so unusually high cells are visible for inspection.
Those flags are not doublet calls; doublet detection should use a
dedicated method.

## Annotation

``` r

annotated_bulk <- annotate_mad_qc(bulk_metadata, bulk_qc)
head(annotated_bulk)
#>   sample_id condition  batch library_size detected_genes mapping_rate
#> 1  bulk_001   control batch1      3000000           4500    0.5200000
#> 2  bulk_002   treated batch2      4285714           5214    0.5542857
#> 3  bulk_003   control batch3      5571429           5929    0.5885714
#> 4  bulk_004   treated batch1      6857143           6643    0.6228571
#> 5  bulk_005   control batch2      8142857           7357    0.6571429
#> 6  bulk_006   treated batch3      9428571           8071    0.6914286
#>   assigned_rate rrna_rate library_size_mad_outlier detected_genes_mad_outlier
#> 1          0.35      0.16                     TRUE                       TRUE
#> 2          0.38      0.18                     TRUE                       TRUE
#> 3          0.41      0.20                     TRUE                       TRUE
#> 4          0.44      0.22                     TRUE                       TRUE
#> 5          0.47      0.24                     TRUE                       TRUE
#> 6          0.50      0.26                     TRUE                       TRUE
#>   mapping_rate_mad_outlier assigned_rate_mad_outlier rrna_rate_mad_outlier
#> 1                     TRUE                      TRUE                  TRUE
#> 2                     TRUE                      TRUE                  TRUE
#> 3                     TRUE                      TRUE                  TRUE
#> 4                     TRUE                      TRUE                  TRUE
#> 5                     TRUE                      TRUE                  TRUE
#> 6                     TRUE                      TRUE                  TRUE
#>   mad_qc_outlier
#> 1           TRUE
#> 2           TRUE
#> 3           TRUE
#> 4           TRUE
#> 5           TRUE
#> 6           TRUE
```

Annotation adds flag columns. It does not filter rows or change the
original QC measurements.

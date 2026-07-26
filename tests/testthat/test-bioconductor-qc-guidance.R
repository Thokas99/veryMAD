test_that("single-cell count metrics are transformed while percent.mt stays raw", {
  cells <- veryMAD:::.simulate_single_cell_qc_metadata()
  metrics <- c(nCount_RNA = "both", nFeature_RNA = "both", percent.mt = "upper")
  transform <- c(nCount_RNA = "log1p", nFeature_RNA = "log1p")

  report <- mad_qc(cells, metrics = metrics, transform = transform)

  count_row <- report[report$metric == "nCount_RNA", ][1, ]
  feature_row <- report[report$metric == "nFeature_RNA", ][1, ]
  mito_row <- report[report$metric == "percent.mt", ][1, ]

  expect_equal(count_row$value, log1p(count_row$raw_value))
  expect_equal(feature_row$value, log1p(feature_row$raw_value))
  expect_equal(mito_row$value, mito_row$raw_value)
  expect_true(any(report$is_outlier[report$metric == "nCount_RNA" & report$direction == "both"]))
  expect_true(any(report$is_outlier[report$metric == "nFeature_RNA" & report$direction == "both"]))
  expect_equal(report$raw_value[report$metric == "nCount_RNA"], cells$nCount_RNA)
})

test_that("bulk count metrics are transformed while rates stay raw", {
  bulk <- veryMAD:::.simulate_bulk_qc_metadata()
  metrics <- c(
    library_size = "lower",
    detected_genes = "lower",
    mapping_rate = "lower",
    assigned_rate = "lower",
    rrna_rate = "upper"
  )
  transform <- c(library_size = "log1p", detected_genes = "log1p")

  report <- mad_qc(bulk, metrics = metrics, transform = transform)

  lib_row <- report[report$metric == "library_size", ][1, ]
  gene_row <- report[report$metric == "detected_genes", ][1, ]
  mapping_row <- report[report$metric == "mapping_rate", ][1, ]
  assigned_row <- report[report$metric == "assigned_rate", ][1, ]
  rrna_row <- report[report$metric == "rrna_rate", ][1, ]

  expect_equal(lib_row$value, log1p(lib_row$raw_value))
  expect_equal(gene_row$value, log1p(gene_row$raw_value))
  expect_equal(mapping_row$value, mapping_row$raw_value)
  expect_equal(assigned_row$value, assigned_row$raw_value)
  expect_equal(rrna_row$value, rrna_row$raw_value)
})

test_that("canonical simulations expose only documented metric columns", {
  bulk <- veryMAD:::.simulate_bulk_qc_metadata()
  cells <- veryMAD:::.simulate_single_cell_qc_metadata()

  expect_equal(
    names(bulk),
    c("sample_id", "condition", "batch", "library_size", "detected_genes",
      "mapping_rate", "assigned_rate", "rrna_rate")
  )
  expect_equal(names(cells), c("cell_id", "nCount_RNA", "nFeature_RNA", "percent.mt"))
  expect_false("doublet_score" %in% names(cells))
  expect_false("duplication_rate" %in% names(bulk))
  expect_false("percent_mitochondrial" %in% names(bulk))
})

test_that("custom metrics default to identity and names do not imply transforms", {
  data <- data.frame(
    library_size = c(10, 20, 30, 40, 100000),
    custom_metric = c(1, 2, 3, 4, 5)
  )

  report <- mad_qc(data, metrics = c(library_size = "both", custom_metric = "both"))

  expect_equal(report$value[report$metric == "library_size"], data$library_size)
  expect_equal(report$value[report$metric == "custom_metric"], data$custom_metric)
})

test_that("Seurat report and annotation preserve explicit flags without filtering", {
  skip_if_not_installed("SeuratObject")

  counts <- matrix(
    rpois(100, lambda = 5),
    nrow = 10,
    dimnames = list(paste0("gene", seq_len(10)), paste0("cell", seq_len(10)))
  )
  object <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = counts))
  object[["percent.mt"]] <- seq(1, 25, length.out = ncol(object))

  metrics <- c(nCount_RNA = "both", nFeature_RNA = "both", percent.mt = "upper")
  transform <- c(nCount_RNA = "log1p", nFeature_RNA = "log1p")

  report <- mad_qc_seurat(object, metrics = metrics, transform = transform, action = "report")
  annotated <- mad_qc_seurat(object, metrics = metrics, transform = transform, action = "annotate")

  expect_setequal(unique(report$metric), names(metrics))
  expect_equal(ncol(annotated), ncol(object))
  expect_equal(SeuratObject::LayerData(annotated, layer = "counts"), SeuratObject::LayerData(object, layer = "counts"))
  expect_true("mad_qc_outlier" %in% names(annotated[[]]))
})

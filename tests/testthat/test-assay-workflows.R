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

single_cell_metrics <- c(
  nCount_RNA = "both",
  nFeature_RNA = "both",
  percent.mt = "upper"
)

single_cell_transform <- c(
  nCount_RNA = "log10p",
  nFeature_RNA = "log10p"
)

test_that("bulk simulation is realistic, related, and reproducible", {
  bulk <- veryMAD:::.simulate_bulk_qc_metadata()
  again <- veryMAD:::.simulate_bulk_qc_metadata()

  expect_equal(nrow(bulk), 240L)
  expect_equal(bulk, again)
  expect_equal(
    names(bulk),
    c("sample_id", "condition", "batch", "library_size", "detected_genes",
      "mapping_rate", "assigned_rate", "rrna_rate")
  )
  expect_true(all(bulk$library_size > 0))
  expect_true(all(bulk$detected_genes > 0))
  expect_true(all(bulk$mapping_rate >= 0 & bulk$mapping_rate <= 1))
  expect_true(all(bulk$assigned_rate >= 0 & bulk$assigned_rate <= 1))
  expect_true(all(bulk$rrna_rate >= 0 & bulk$rrna_rate <= 1))
})

test_that("bulk workflow flags deliberately poor libraries", {
  bulk <- veryMAD:::.simulate_bulk_qc_metadata()
  report <- mad_qc(bulk, bulk_metrics, transform = bulk_transform)
  summary <- summarize_mad_qc(report, level = "observation")

  expect_true(any(summary$mad_qc_outlier[seq_len(8L)]))
  expect_true(any(report$is_outlier[report$metric == "library_size"]))
  expect_true(any(report$is_outlier[report$metric == "detected_genes"]))
  expect_true(any(report$is_outlier[report$metric == "rrna_rate"]))
})

test_that("single-cell simulation is reproducible and omits doublet scores", {
  cells <- veryMAD:::.simulate_single_cell_qc_metadata()
  again <- veryMAD:::.simulate_single_cell_qc_metadata()

  expect_equal(nrow(cells), 1200L)
  expect_equal(cells, again)
  expect_equal(names(cells), c("cell_id", "nCount_RNA", "nFeature_RNA", "percent.mt"))
  expect_false("doublet_score" %in% names(cells))
  expect_true(all(cells$nCount_RNA > 0))
  expect_true(all(cells$nFeature_RNA > 0))
  expect_true(all(cells$percent.mt >= 0))
}
)

test_that("single-cell workflow uses explicit log count transforms and raw percent.mt", {
  cells <- veryMAD:::.simulate_single_cell_qc_metadata()
  report <- mad_qc(cells, single_cell_metrics, transform = single_cell_transform)

  expect_true(any(report$is_outlier[report$metric == "nCount_RNA"]))
  expect_true(any(report$is_outlier[report$metric == "nFeature_RNA"]))
  expect_true(any(report$is_outlier[report$metric == "percent.mt"]))
  expect_equal(report$value[report$metric == "percent.mt"], report$raw_value[report$metric == "percent.mt"])
}
)

test_that("grouped bulk QC keeps rates on the raw scale", {
  bulk <- veryMAD:::.simulate_bulk_qc_metadata()
  grouped <- mad_qc(bulk, c(mapping_rate = "lower"), group_by = "batch")

  expect_equal(grouped$value, grouped$raw_value)
  expect_true(all(grouped$metric == "mapping_rate"))
}
)

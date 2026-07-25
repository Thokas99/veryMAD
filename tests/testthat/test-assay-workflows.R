bulk_metrics <- c(library_size = "lower", detected_genes = "lower",
  mapping_rate = "lower", duplication_rate = "upper",
  percent_mitochondrial = "upper")

cell_metrics <- c(nCount_RNA = "lower", nFeature_RNA = "lower",
  percent.mt = "upper", doublet_score = "upper")

test_that("bulk simulation is realistic, related, and reproducible", {
  bulk <- veryMAD:::.simulate_bulk_qc_metadata()
  expect_equal(nrow(bulk), 150L)
  expect_named(bulk, c("sample_id", "condition", "batch", "library_size",
    "detected_genes", "mapping_rate", "duplication_rate", "percent_mitochondrial"))
  expect_equal(length(unique(bulk$sample_id)), nrow(bulk))
  expect_equal(sort(unique(bulk$condition)), c("Control", "Treatment"))
  expect_equal(sort(unique(bulk$batch)), c("Batch1", "Batch2", "Batch3"))
  complete <- complete.cases(bulk[c("library_size", "detected_genes")])
  expect_gt(stats::cor(bulk$library_size[complete], bulk$detected_genes[complete]), 0.3)
  expect_true(all(bulk$mapping_rate[!is.na(bulk$mapping_rate)] >= 0 & bulk$mapping_rate[!is.na(bulk$mapping_rate)] <= 1))
  expect_true(all(bulk$duplication_rate[!is.na(bulk$duplication_rate)] >= 0 & bulk$duplication_rate[!is.na(bulk$duplication_rate)] <= 1))
  expect_true(all(bulk$percent_mitochondrial[!is.na(bulk$percent_mitochondrial)] >= 0 & bulk$percent_mitochondrial[!is.na(bulk$percent_mitochondrial)] <= 100))
  expect_true(anyNA(bulk))
  expect_equal(bulk, veryMAD:::.simulate_bulk_qc_metadata())
})

test_that("bulk QC is globally thresholded and every metric has failures", {
  bulk <- veryMAD:::.simulate_bulk_qc_metadata()
  report <- mad_qc(bulk, bulk_metrics,
    transform = c(library_size = "log10p", detected_genes = "log10p"))
  expect_null(attr(report, "group_by"))
  expect_equal(report$id[report$metric == "library_size"], bulk$sample_id)
  expect_false(any(c("condition", "batch") %in% names(report)))
  expect_true(all(vapply(split(report$median, report$metric), function(x) length(unique(x)) == 1L, logical(1))))
  expect_true(all(vapply(split(report$is_outlier, report$metric), function(x) any(x %in% TRUE), logical(1))))
  expect_true(any(report$is_outlier %in% FALSE)); expect_true(any(report$is_outlier %in% TRUE)); expect_true(anyNA(report$is_outlier))
})

test_that("single-cell simulation respects ranges and relationships", {
  cells <- veryMAD:::.simulate_single_cell_qc_metadata()
  expect_equal(nrow(cells), 1000L)
  expect_named(cells, c("cell_id", "nCount_RNA", "nFeature_RNA", "percent.mt", "doublet_score"))
  expect_equal(length(unique(cells$cell_id)), nrow(cells))
  complete <- complete.cases(cells[c("nCount_RNA", "nFeature_RNA")])
  expect_gt(stats::cor(cells$nCount_RNA[complete], cells$nFeature_RNA[complete]), 0.3)
  expect_true(all(cells$nCount_RNA[!is.na(cells$nCount_RNA)] >= 0))
  expect_true(all(cells$nFeature_RNA[complete] <= cells$nCount_RNA[complete]))
  expect_true(all(cells$percent.mt[!is.na(cells$percent.mt)] >= 0 & cells$percent.mt[!is.na(cells$percent.mt)] <= 100))
  expect_true(all(cells$doublet_score[!is.na(cells$doublet_score)] >= 0 & cells$doublet_score[!is.na(cells$doublet_score)] <= 1))
  expect_true(anyNA(cells))
})

test_that("single-cell QC is ungrouped and every metric has failures", {
  cells <- veryMAD:::.simulate_single_cell_qc_metadata()
  report <- mad_qc(cells, cell_metrics,
    transform = c(nCount_RNA = "log10p", nFeature_RNA = "log10p"))
  expect_null(attr(report, "group_by"))
  expect_false(any(c("sample_id", "orig.ident") %in% names(cells)))
  expect_true(all(vapply(split(report$is_outlier, report$metric), function(x) any(x %in% TRUE), logical(1))))
  expect_true(any(report$is_outlier %in% FALSE)); expect_true(any(report$is_outlier %in% TRUE)); expect_true(anyNA(report$is_outlier))
})

test_that("grouped thresholds exist only when explicitly requested", {
  bulk <- veryMAD:::.simulate_bulk_qc_metadata()
  ungrouped <- mad_qc(bulk, c(mapping_rate = "lower"))
  grouped <- mad_qc(bulk, c(mapping_rate = "lower"), group_by = "batch")
  expect_null(attr(ungrouped, "group_by"))
  expect_equal(attr(grouped, "group_by"), "batch")
  expect_equal(length(unique(ungrouped$median)), 1L)
  expect_gt(length(unique(grouped$median)), 1L)
})

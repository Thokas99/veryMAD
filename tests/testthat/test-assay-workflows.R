test_that("bulk metadata uses one global reference distribution by default", {
  bulk <- veryMAD:::.simulate_bulk_qc_metadata()
  expect_equal(nrow(bulk), 250L)
  report <- mad_qc(
    bulk,
    metrics = c(library_size = "lower", detected_genes = "lower",
      mapping_rate = "lower", duplication_rate = "upper",
      percent_mitochondrial = "upper"),
    transform = c(library_size = "log10p", detected_genes = "log10p")
  )
  expect_null(attr(report, "group_by"))
  expect_equal(report$id[report$metric == "library_size"], bulk$sample_id)
  expect_false("condition" %in% names(report))
  expect_true("condition" %in% names(attr(report, "observation_metadata")))
  expect_equal(length(unique(report$median[report$metric == "library_size"])), 1L)
  expect_equal(length(unique(report$mad[report$metric == "library_size"])), 1L)
  expect_equal(length(unique(report$lower[report$metric == "library_size"])), 1L)
  flags <- report$is_outlier
  expect_true(any(flags %in% FALSE))
  expect_true(any(flags %in% TRUE))
  expect_true(anyNA(flags))
})

test_that("grouped thresholds exist only when explicitly requested", {
  bulk <- veryMAD:::.simulate_bulk_qc_metadata()
  ungrouped <- mad_qc(bulk, c(mapping_rate = "lower"))
  grouped <- mad_qc(bulk, c(mapping_rate = "lower"), group_by = "condition")
  expect_null(attr(ungrouped, "group_by"))
  expect_equal(attr(grouped, "group_by"), "condition")
  expect_equal(length(unique(ungrouped$median)), 1L)
  expect_gt(length(unique(grouped$median)), 1L)
})

test_that("single-cell metadata remains ungrouped and includes all QC states", {
  cells <- veryMAD:::.simulate_single_cell_qc_metadata()
  expect_equal(nrow(cells), 600L)
  expect_false(any(c("sample_id", "orig.ident") %in% names(cells)))
  report <- mad_qc(cells, c(nCount_RNA = "lower", nFeature_RNA = "lower",
    percent.mt = "upper", doublet_score = "upper"),
    transform = c(nCount_RNA = "log10p", nFeature_RNA = "log10p"))
  expect_null(attr(report, "group_by"))
  expect_equal(report$id[report$metric == "nCount_RNA"], cells$cell_id)
  expect_true(any(report$is_outlier %in% FALSE))
  expect_true(any(report$is_outlier %in% TRUE))
  expect_true(anyNA(report$is_outlier))
})

test_that("calculation grouping and visual faceting stay independent", {
  skip_if_not_installed("ggplot2")
  bulk <- veryMAD:::.simulate_bulk_qc_metadata()
  global <- mad_qc(bulk, c(mapping_rate = "lower"))
  global_plot <- plot_mad_qc(global)
  faceted_plot <- plot_mad_qc(global, facet_by = "condition")
  expect_equal(nlevels(global_plot$data$threshold_group), 1L)
  expect_equal(nlevels(global_plot$data$plot_facet), 1L)
  expect_equal(nlevels(faceted_plot$data$threshold_group), 1L)
  expect_equal(nlevels(faceted_plot$data$plot_facet), 2L)

  grouped <- mad_qc(bulk, c(mapping_rate = "lower"), group_by = "condition")
  grouped_plot <- plot_mad_qc(grouped)
  expect_equal(nlevels(grouped_plot$data$threshold_group), 2L)
  expect_equal(nlevels(grouped_plot$data$plot_facet), 1L)
  grouped_thresholds <- veryMAD:::.qc_plot_thresholds(grouped_plot$data, FALSE)
  expect_gt(length(unique(grouped_thresholds$threshold)), 1L)
})

plot_fixture <- function() {
  d <- data.frame(sample = rep(c("A", "B"), each = 4),
                  x = c(1, 2, 3, 20, 2, 3, 4, 30),
                  y = c(2, 3, 20, 4, 2, 30, 4, 5))
  mad_qc(d, c(x = "upper", y = "both"), nmads = 1, group_by = "sample")
}

test_that("plot_mad_qc returns plots for both views and metric selection", {
  skip_if_not_installed("ggplot2")
  qc <- plot_fixture(); original <- qc
  distribution <- plot_mad_qc(qc, metrics = "x")
  index <- plot_mad_qc(qc, type = "index", facet_by = "sample")
  expect_s3_class(distribution, "ggplot")
  expect_s3_class(index, "ggplot")
  expect_equal(unique(as.character(distribution$data$metric)), "x")
  expect_equal(nlevels(distribution$data$threshold_group), 2L)
  expect_equal(qc, original)
  expect_error(plot_mad_qc(qc, metrics = "missing"), "Unknown metric")
})

test_that("plot thresholds omit infinite and missing limits", {
  skip_if_not_installed("ggplot2")
  qc <- plot_fixture()
  qc$lower[1] <- NA_real_
  qc$upper[2] <- Inf
  qc$threshold_group <- veryMAD:::.qc_plot_labels(qc, attr(qc, "group_by"), "All observations")
  qc$plot_facet <- veryMAD:::.qc_plot_labels(qc, NULL, "All observations")
  thresholds <- veryMAD:::.qc_plot_thresholds(qc, by_observation = TRUE)
  expect_true(all(is.finite(thresholds$threshold)))
  expect_s3_class(plot_mad_qc(qc, show_thresholds = FALSE, show_legend = FALSE), "ggplot")
})

test_that("empty plots and optional dependency failure are deliberate", {
  skip_if_not_installed("ggplot2")
  empty <- mad_qc(data.frame(x = numeric()), c(x = "both"))
  expect_s3_class(plot_mad_qc(empty), "ggplot")
  testthat::local_mocked_bindings(
    .require_namespace = function(package, feature) stop("ggplot2 unavailable", call. = FALSE),
    .package = "veryMAD"
  )
  expect_error(plot_mad_qc(plot_fixture()), "ggplot2 unavailable")
})

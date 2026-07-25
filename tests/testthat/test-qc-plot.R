plot_fixture <- function(grouped = TRUE) {
  d <- data.frame(sample = rep(c("A", "B"), each = 5),
    x = c(1, 2, 3, 20, NA, 2, 3, 4, 30, 5),
    y = c(2, 3, 20, 4, 5, 2, 30, 4, 5, 6))
  mad_qc(d, c(x = "upper", y = "both"), nmads = 1,
    group_by = if (grouped) "sample" else NULL)
}

test_that("plot_mad_qc returns both views without modifying input", {
  skip_if_not_installed("ggplot2")
  qc <- plot_fixture(); original <- qc
  distribution <- plot_mad_qc(qc, metrics = "x")
  index <- plot_mad_qc(qc, type = "index", facet_by = "sample")
  expect_s3_class(distribution, "ggplot")
  expect_s3_class(index, "ggplot")
  expect_equal(unique(as.character(distribution$data$metric)), "x")
  expect_equal(qc, original)
  expect_error(plot_mad_qc(qc, metrics = "missing"), "Unknown metric")
})

test_that("plot status colors and missing indicators are explicit", {
  skip_if_not_installed("ggplot2")
  qc <- plot_fixture()
  for (type in c("distribution", "index")) {
    plot <- plot_mad_qc(qc, type = type)
    colour_scale <- plot$scales$get_scales("colour")
    expect_equal(unname(colour_scale$palette(3)), c("#009E73", "#0072B2", "#D55E00"))
    expect_equal(levels(plot$data$status), c("Pass", "Not evaluated", "Outlier"))
    expect_match(plot$labels$caption, "Blue triangle")
    built <- expect_warning(ggplot2::ggplot_build(plot), NA)
    expect_true(any(vapply(built$data, function(layer) any(layer$shape %in% 17), logical(1))))
  }
})

test_that("thresholds are active, styled, and separated by group", {
  skip_if_not_installed("ggplot2")
  qc <- plot_fixture()
  plot <- plot_mad_qc(qc)
  thresholds <- veryMAD:::.qc_plot_thresholds(plot$data, by_observation = TRUE)
  expect_true(all(is.finite(thresholds$threshold)))
  expect_false(any(thresholds$limit[thresholds$metric == "x"] == "Lower MAD threshold"))
  expect_equal(levels(thresholds$limit), c("Lower MAD threshold", "Upper MAD threshold"))
  expect_gt(length(unique(thresholds$trajectory)), 2L)
  linetype_scale <- plot$scales$get_scales("linetype")
  expect_equal(unname(linetype_scale$palette(2)), c("dashed", "dotdash"))
  expect_s3_class(plot_mad_qc(qc, show_thresholds = FALSE, show_legend = FALSE), "ggplot")
  expect_equal(plot_mad_qc(qc, show_legend = FALSE)$theme$legend.position, "none")
})

test_that("calculation grouping and visual faceting remain independent", {
  skip_if_not_installed("ggplot2")
  global <- plot_mad_qc(plot_fixture(grouped = FALSE))
  faceted <- plot_mad_qc(plot_fixture(grouped = FALSE), facet_by = "sample")
  grouped <- plot_mad_qc(plot_fixture(grouped = TRUE))
  expect_equal(nlevels(global$data$threshold_group), 1L)
  expect_equal(nlevels(faceted$data$threshold_group), 1L)
  expect_equal(nlevels(faceted$data$plot_facet), 2L)
  expect_equal(nlevels(grouped$data$threshold_group), 2L)
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

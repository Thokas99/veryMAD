test_that("batch-specific thresholds apply to their own observations", {
  d <- data.frame(
    value = c(1:5, 101:105, 1000),
    sample = factor(c(rep("A", 5), rep("B", 6)), levels = c("A", "B", "unused"))
  )
  report <- mad_qc(d, c(value = "upper"), batch = d$sample,
                   output = "report", verbose = FALSE)

  expect_equal(report$thresholds$batch, c("A", "B"))
  expect_equal(report$thresholds$median, c(3, 103.5))
  expect_equal(report$thresholds$n_used, c(5L, 6L))
  expect_equal(report$thresholds$n_flagged, c(0L, 1L))
  expect_true(tail(report$flags$value, 1) %in% TRUE)
})

test_that("reference rows estimate thresholds while all rows are flagged", {
  d <- data.frame(value = c(1:5, 1000))
  report <- mad_qc(d, c(value = "upper"), subset = c(rep(TRUE, 5), FALSE),
                   output = "report", verbose = FALSE)

  expect_equal(report$thresholds$median, 3)
  expect_equal(report$thresholds$n_used, 5L)
  expect_equal(report$settings$n_reference, 5L)
  expect_true(report$flags$value[6] %in% TRUE)
})

test_that("undersized batches retain unresolved flags", {
  d <- data.frame(value = c(1:5, 20, 21), sample = c(rep("A", 5), rep("B", 2)))
  expect_warning(report <- mad_qc(
    d, c(value = "upper"), batch = d$sample, min_n = 3,
    output = "report", verbose = FALSE
  ), "value \\[B\\]")

  expect_equal(report$thresholds$status, c("ok", "insufficient_n"))
  expect_true(all(is.na(report$flags$value[6:7])))
  expect_equal(sum(report$thresholds$status != "ok"), 1L)
  expect_output(print(report), "Metrics without estimated thresholds: 1")
})

test_that("batch and reference selectors are validated", {
  d <- data.frame(value = 1:5)
  expect_error(mad_qc(d, c(value = "upper"), batch = "one", verbose = FALSE), "batch")
  expect_error(mad_qc(d, c(value = "upper"), batch = c("A", NA, "A", "A", "A"), verbose = FALSE), "batch")
  expect_error(mad_qc(d, c(value = "upper"), subset = c(TRUE, NA, TRUE, TRUE, TRUE), verbose = FALSE), "subset")
  expect_error(mad_qc(d, c(value = "upper"), subset = 1:5, verbose = FALSE), "subset")
})

test_that("plot_mad_qc shows QC status by batch", {
  skip_if_not_installed("ggplot2")
  d <- data.frame(value = c(1:5, 100), sample = rep(c("A", "B"), each = 3))
  annotated <- mad_qc(d, c(value = "upper"), batch = d$sample,
                      min_n = 3, verbose = FALSE)
  plot <- plot_mad_qc(annotated, "value", batch = d$sample)

  expect_s3_class(plot, "ggplot")
  expect_equal(plot$data$batch, factor(d$sample, levels = c("A", "B")))
  expect_true(any(plot$data$status == "Flagged"))
  expect_true(any(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomViolin"), logical(1))))
  expect_true(any(vapply(plot$layers, function(layer) inherits(layer$geom, "GeomBoxplot"), logical(1))))
})

test_that("plot_mad_qc handles single-observation batches without plot warnings", {
  skip_if_not_installed("ggplot2")
  d <- data.frame(value = 1:6, sample = c(rep("A", 5), "B"))
  expect_warning(
    annotated <- mad_qc(d, c(value = "upper"), batch = d$sample, verbose = FALSE),
    "Insufficient usable observations"
  )
  plot <- plot_mad_qc(annotated, "value", batch = d$sample)

  expect_no_warning(ggplot2::ggplotGrob(plot))
})

test_that("plot_mad_qc_scatter reports within-batch Spearman rho and n", {
  skip_if_not_installed("ggplot2")
  d <- data.frame(
    x = c(1, 2, NA, 4, 5, 5, 5),
    y = c(4, 3, 2, 1, 3, 1, 2),
    sample = c(rep("A", 4), rep("B", 3)),
    x_mad_outlier = c(FALSE, TRUE, NA, FALSE, FALSE, FALSE, FALSE),
    y_mad_outlier = rep(FALSE, 7)
  )

  plot <- plot_mad_qc_scatter(d, "x", "y", batch = d$sample)
  labels <- plot$layers[[2]]$data

  expect_s3_class(plot, "ggplot")
  expect_equal(labels$rho[1], -1)
  expect_true(is.na(labels$rho[2]))
  expect_equal(labels$n, c(3L, 3L))
  expect_equal(
    ggplot2::ggplot_build(plot)$data[[2]]$label,
    c("rho = -1.00\nn = 3", "rho = NA\nn = 3")
  )
  expect_no_warning(ggplot2::ggplotGrob(plot))
  expect_true(any(plot$data$status == "Flagged"))
  expect_true(any(plot$data$status == "Unresolved"))
  expect_error(plot_mad_qc_scatter(d, "x", "x"), "distinct")
})

test_that("plot_mad_qc_scatter overlays raw report thresholds by batch", {
  skip_if_not_installed("ggplot2")
  d <- data.frame(
    x = c(1:5, 101:105),
    y = c(5:1, 110:106),
    sample = rep(c("A", "B"), each = 5)
  )
  metrics <- c(x = "lower", y = "upper")
  transforms <- c(x = "log1p", y = "log1p")
  annotated <- mad_qc(d, metrics, transform = transforms,
                      batch = d$sample, verbose = FALSE)
  report <- mad_qc(d, metrics, transform = transforms, batch = d$sample,
                   output = "report", verbose = FALSE)

  plot <- plot_mad_qc_scatter(
    annotated, "x", "y", batch = d$sample, report = report
  )
  vline <- plot$layers[[which(vapply(plot$layers, function(layer) {
    inherits(layer$geom, "GeomVline")
  }, logical(1)))]]
  hline <- plot$layers[[which(vapply(plot$layers, function(layer) {
    inherits(layer$geom, "GeomHline")
  }, logical(1)))]]

  expect_equal(vline$data$cutoff, report$thresholds$lower_raw[report$thresholds$metric == "x"])
  expect_equal(hline$data$cutoff, report$thresholds$upper_raw[report$thresholds$metric == "y"])
  expect_equal(as.character(vline$data$batch), report$thresholds$batch[report$thresholds$metric == "x"])
  expect_no_warning(ggplot2::ggplotGrob(plot))
})

test_that("plot_mad_qc_scatter omits unresolved limits and repeats pooled limits", {
  skip_if_not_installed("ggplot2")
  d <- data.frame(x = rep(1, 5), y = 1:5,
                  sample = rep(c("A", "B"), length.out = 5))
  metrics <- c(x = "lower", y = "upper")
  annotated <- mad_qc(d, metrics, verbose = FALSE)
  pooled_report <- mad_qc(d, metrics, output = "report", verbose = FALSE)

  pooled_plot <- plot_mad_qc_scatter(
    annotated, "x", "y", batch = d$sample, report = pooled_report
  )
  geoms <- vapply(pooled_plot$layers, function(layer) class(layer$geom)[1], character(1))
  expect_false("GeomVline" %in% geoms)
  expect_true("GeomHline" %in% geoms)
  hline <- pooled_plot$layers[[which(geoms == "GeomHline")]]
  expect_equal(nrow(hline$data), length(unique(d$sample)))
  expect_true(all(hline$data$cutoff == pooled_report$thresholds$upper_raw[
    pooled_report$thresholds$metric == "y"
  ]))
})

test_that("plot_mad_qc_scatter requires a matching QC report", {
  skip_if_not_installed("ggplot2")
  d <- data.frame(x = c(1:5, 100), y = 6:1)
  annotated <- mad_qc(d, c(x = "lower", y = "upper"), verbose = FALSE)
  report <- mad_qc(d[5:1, ], c(x = "lower", y = "upper"),
                   output = "report", verbose = FALSE)
  mismatched_report <- mad_qc(d, c(x = "both", y = "upper"),
                              output = "report", verbose = FALSE)

  expect_error(plot_mad_qc_scatter(annotated, "x", "y", report = list()),
               "`verymad_qc` report")
  expect_error(plot_mad_qc_scatter(annotated, "x", "y", report = report),
               "same observations and order")
  expect_error(plot_mad_qc_scatter(annotated, "x", "y", report = mismatched_report),
               "flags must match")
})

test_that("QC plotting rejects empty annotated metadata", {
  skip_if_not_installed("ggplot2")
  d <- data.frame(
    x = numeric(), y = numeric(),
    x_mad_outlier = logical(), y_mad_outlier = logical()
  )

  expect_error(plot_mad_qc(d, "x"), "at least one observation")
  expect_error(plot_mad_qc_scatter(d, "x", "y"), "at least one observation")
})

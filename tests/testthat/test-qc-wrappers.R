wrapper_data <- function() data.frame(
  batch = c("a", "a", "a", "a", "b", "b", "b", "b"),
  count = c(NA, 10, 11, 100, 20, 21, 22, 200),
  rate = c(.9, .91, .92, .1, .8, .81, .82, .2),
  label = letters[1:8], row.names = paste0("obs", 1:8)
)

test_that("wrappers require explicit valid metrics", {
  d <- wrapper_data()
  for (bad in list(NULL, character(), c("lower"), c(x = "lower", x = "upper"))) {
    expect_error(mad_qc_bulk(d, bad, verbose = FALSE), "metrics")
  }
  expect_error(mad_qc_sc(d, verbose = FALSE), "will not guess")
  expect_error(mad_qc_bulk(d, c(missing = "lower"), verbose = FALSE), "Missing")
  expect_error(mad_qc_bulk(d, c(label = "lower"), verbose = FALSE), "numeric")
  expect_error(mad_qc_bulk(d, c(count = "sideways"), verbose = FALSE), "directions")
  expect_identical(names(mad_qc_bulk(d, c(rate = "lower"), verbose = FALSE)),
    c(names(d), "rate_mad_outlier"))
})

test_that("transform defaults and partial overrides are explicit", {
  d <- wrapper_data(); metrics <- c(count = "upper", rate = "lower")
  default <- mad_qc_sc(d, metrics, action = "report", verbose = FALSE)
  expect_equal(default$settings$transform, c(count = "log1p", rate = "log1p"))
  partial <- mad_qc_sc(d, metrics, transform = c(rate = "none"), action = "report", verbose = FALSE)
  expect_equal(partial$settings$transform, c(count = "log1p", rate = "none"))
  none <- mad_qc_sc(d, metrics, transform = "none", action = "report", verbose = FALSE)
  identity <- mad_qc_sc(d, metrics, transform = "identity", action = "report", verbose = FALSE)
  expect_equal(none$thresholds, identity$thresholds)
  expect_equal(mad_qc_sc(transform(d, count = count + 1), metrics, transform = "log10",
    action = "report", verbose = FALSE)$settings$transform, c(count = "log10", rate = "log10"))
  expect_error(mad_qc_sc(transform(d, count = replace(count, 2, -1)), metrics,
    action = "report", verbose = FALSE), "nonnegative")
  expect_error(mad_qc_sc(transform(d, count = replace(count, 2, Inf)), metrics,
    action = "report", verbose = FALSE), "non-finite")
  expect_error(mad_qc_sc(d, metrics, transform = c(other = "none"),
    action = "report", verbose = FALSE), "selected metrics")
  row <- partial$thresholds[partial$thresholds$metric == "count", ]
  expect_equal(row$upper_raw, expm1(row$upper))
})

test_that("bulk handles supported inputs, groups, identifiers, and collisions", {
  d <- wrapper_data(); metrics <- c(count = "upper", rate = "lower")
  for (x in list(d, tibble::as_tibble(d), as.matrix(d[c("count", "rate")]))) {
    out <- mad_qc_bulk(x, metrics, transform = "none", verbose = FALSE)
    expect_true(all(c("count_mad_outlier", "rate_mad_outlier") %in% names(out)))
    expect_equal(nrow(out), 8)
  }
  out <- mad_qc_bulk(d, metrics, group_by = "batch", transform = "none", verbose = FALSE)
  report <- attr(out, "mad_qc")
  expect_s3_class(report, "mad_qc_result")
  expect_equal(nrow(report$thresholds), 4)
  expect_equal(report$settings$observation_ids, rownames(d))
  expect_equal(rownames(out), rownames(d))
  expect_error(mad_qc_bulk(out, metrics, verbose = FALSE), "already exist")
  expect_silent(mad_qc_bulk(out, metrics, verbose = FALSE, overwrite = TRUE))
})

test_that("empty and zero-MAD inputs have neutral tested behavior", {
  empty <- mad_qc_bulk(data.frame(x = numeric()), c(x = "both"), verbose = FALSE)
  expect_equal(nrow(empty), 0); expect_equal(nrow(attr(empty, "mad_qc")$thresholds), 0)
  report <- mad_qc_sc(data.frame(x = c(1, 1, 1)), c(x = "both"),
    transform = "none", action = "report", verbose = FALSE)
  expect_false(any(report$flags$x, na.rm = TRUE))
  expect_equal(report$thresholds$mad, 0)
})

test_that("single-cell annotation exposes one three-valued flag and never filters", {
  d <- wrapper_data(); metrics <- c(count = "upper", rate = "lower")
  before <- d
  report <- mad_qc_sc(d, metrics, transform = "none", action = "report", verbose = FALSE)
  expect_identical(d, before)
  out <- mad_qc_sc(d, metrics, transform = "none", verbose = FALSE)
  expect_identical(setdiff(names(out), names(d)), "mad_qc_outlier")
  expect_equal(nrow(out), nrow(d)); expect_true(is.na(out$mad_qc_outlier[1]))
  expect_error(mad_qc_sc(out, metrics, verbose = FALSE), "already exists")
  expect_silent(mad_qc_sc(out, metrics, verbose = FALSE, overwrite = TRUE))
  expect_equal(out$mad_qc_outlier, veryMAD:::.overall_qc_flag(report$flags[names(metrics)]))
})

test_that("grades use deterministic unique-observation proportions", {
  expect_equal(veryMAD:::.qc_grade(0), "Wonderful! :D")
  expect_match(veryMAD:::.qc_grade(.05), "Looking good", fixed = TRUE)
  expect_equal(veryMAD:::.qc_grade(.05001), "Not so good...")
  expect_equal(veryMAD:::.qc_grade(.15), "Not so good...")
  expect_match(veryMAD:::.qc_grade(.15001), "wet lab", fixed = TRUE)
})

test_that("CLI summaries are concise, raw-scale, and suppressible", {
  d <- wrapper_data(); metrics <- c(count = "upper", rate = "lower")
  expect_silent(mad_qc_bulk(d, metrics, verbose = FALSE))
  text <- capture_messages(mad_qc_bulk(d, metrics, transform = "none", verbose = TRUE))
  text <- paste(text, collapse = "\n")
  expect_match(text, "8 observations checked across 2 explicit metrics")
  expect_match(text, "flagged")
  grouped <- capture_messages(mad_qc_bulk(d, metrics, group_by = "batch", transform = "none", verbose = TRUE))
  expect_lt(length(grouped), 20)
  many <- data.frame(group = seq_len(12), x = seq_len(12))
  many_text <- capture_messages(mad_qc_bulk(many, c(x = "both"), group_by = "group", verbose = TRUE))
  expect_lt(length(many_text), 20)
  expect_true(any(grepl("threshold rows", many_text)))
})

test_that("startup banner is plain, dynamic, and quiet option is respected", {
  banner <- veryMAD:::.verymad_banner("9.8.7")
  expect_match(banner, "__  __")
  expect_match(banner, "veryMAD 9.8.7")
  expect_false(grepl("animal|face|lab drawing", banner, ignore.case = TRUE))
  old <- options(veryMAD.quiet = TRUE); on.exit(options(old), add = TRUE)
  expect_true(isTRUE(getOption("veryMAD.quiet")))
  expect_silent(veryMAD:::.verymad_startup_message())
  expect_message(veryMAD:::.verymad_startup_message(FALSE, "9.8.7"), "veryMAD 9.8.7")
  expect_silent(suppressPackageStartupMessages(veryMAD:::.verymad_startup_message(FALSE, "9.8.7")))
})

test_that("Seurat annotation preserves object and adds only combined flag", {
  skip_if_not_installed("SeuratObject")
  object <- suppressWarnings(SeuratObject::CreateSeuratObject(matrix(seq_len(32), nrow = 4,
    dimnames = list(paste0("g", 1:4), paste0("c", 1:8)))))
  object$rate <- wrapper_data()$rate
  before_cells <- colnames(object); before_assay <- SeuratObject::DefaultAssay(object)
  out <- mad_qc_sc(object, c(nCount_RNA = "lower", rate = "lower"),
    transform = c(rate = "none"), verbose = FALSE)
  expect_equal(colnames(out), before_cells)
  expect_equal(SeuratObject::DefaultAssay(out), before_assay)
  expect_true("mad_qc_outlier" %in% names(out[[]]))
  expect_false("rate_mad_outlier" %in% names(out[[]]))
})

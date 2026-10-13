test_that("mad_qc annotates explicit tails without filtering or changing values", {
  d <- qc_data()
  before <- d
  out <- mad_qc(d, c(lower = "lower", upper = "upper", both = "both"), verbose = FALSE)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), nrow(d))
  expect_identical(rownames(out), rownames(d))
  expect_identical(out[ names(d) ], d)
  expect_true(all(c("lower_mad_outlier", "upper_mad_outlier", "both_mad_outlier", "mad_qc_outlier") %in% names(out)))
  expect_true(any(out$lower_mad_outlier))
  expect_true(any(out$upper_mad_outlier))
  expect_true(any(out$both_mad_outlier))
  expect_true(any(out$mad_qc_outlier))
  expect_identical(d, before)
})

test_that("report output has the stable compact contract", {
  r <- mad_qc(qc_data(), c(lower = "lower", upper = "upper"), output = "report", verbose = FALSE)
  expect_s3_class(r, "verymad_qc")
  expect_identical(names(r), c("flags", "thresholds", "settings"))
  expect_identical(names(r$flags), c("id", "lower", "upper", "mad_qc_outlier"))
  expect_identical(names(r$thresholds), c("metric", "direction", "transform", "median", "mad", "lower", "upper", "lower_raw", "upper_raw", "status", "n_used", "n_flagged", "n_unresolved"))
  expect_true(all(vapply(r$flags[-1], is.logical, logical(1))))
  expect_identical(names(r$settings), c("metrics", "transform", "nmads", "constant", "min_n", "zero_mad"))
  expect_false("verbose" %in% names(r$settings))
  expect_false("annotated" %in% names(r))
})

test_that("combined flags use three-valued any logic and are always present", {
  d <- data.frame(a = c(1:4, 100), missing = rep(NA_real_, 5))
  r <- mad_qc(d, c(a = "upper", missing = "upper"), output = "report", min_n = 3, verbose = FALSE)
  expect_true(any(r$flags$mad_qc_outlier %in% TRUE))
  expect_true(anyNA(r$flags$mad_qc_outlier))

  all_false <- mad_qc(data.frame(a = 1:5), c(a = "both"), output = "report", min_n = 3, verbose = FALSE)
  expect_true(all_false$flags$mad_qc_outlier[1] %in% FALSE)
  expect_false(anyNA(all_false$flags$mad_qc_outlier))

  empty <- mad_qc(data.frame(a = numeric()), c(a = "upper"), output = "report", verbose = FALSE)
  expect_identical(empty$flags$mad_qc_outlier, logical())
})

test_that("combined flags match three-valued any logic for all logical combinations", {
  reference <- function(flags) {
    vapply(seq_len(nrow(flags)), function(i) {
      x <- unlist(flags[i, , drop = FALSE], use.names = FALSE)
      if (any(x %in% TRUE)) TRUE else if (anyNA(x)) NA else FALSE
    }, logical(1))
  }
  values <- expand.grid(a = c(FALSE, TRUE, NA), b = c(FALSE, TRUE, NA),
                        c = c(FALSE, TRUE, NA))

  expect_identical(veryMAD:::.combine_flags(values), reference(values))
  expect_identical(veryMAD:::.combine_flags(values[FALSE, ]), logical())
  expect_identical(veryMAD:::.combine_flags(data.frame(row.names = 1:3)),
                   reference(data.frame(row.names = 1:3)))
})

test_that("threshold reports count finite, flagged, and unresolved observations", {
  d <- data.frame(
    usable = c(1, 2, 3, 4, 100, NA_real_),
    insufficient = c(1, 2, NA_real_, NA_real_, NA_real_, NA_real_),
    missing = rep(NA_real_, 6)
  )
  expect_warning(r <- mad_qc(
    d, c(usable = "upper", insufficient = "lower", missing = "upper"),
    output = "report", min_n = 5, verbose = FALSE
  ), "Insufficient")

  counts <- r$thresholds[, c("metric", "status", "n_used", "n_flagged", "n_unresolved")]
  expect_equal(counts$status, c("ok", "insufficient_n", "all_missing"))
  expect_equal(counts$n_used, c(5L, 2L, 0L))
  expect_equal(counts$n_flagged, c(1L, 0L, 0L))
  expect_equal(counts$n_unresolved, c(1L, 6L, 6L))
  expect_true(r$flags$mad_qc_outlier[5])
  expect_true(is.na(r$flags$mad_qc_outlier[1]))
})

test_that("report identifiers preserve valid observation row names", {
  d <- qc_data()
  rownames(d) <- paste0("cell_", seq_len(nrow(d)))
  report <- mad_qc(d, c(lower = "lower"), output = "report", verbose = FALSE)
  annotated <- mad_qc(d, c(lower = "lower"), verbose = FALSE)

  expect_identical(report$flags$id, rownames(d))
  expect_identical(rownames(annotated), rownames(d))
  expect_identical(annotated$lower, d$lower)
})

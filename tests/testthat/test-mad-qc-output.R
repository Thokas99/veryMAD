test_that("mad_qc annotates explicit tails without filtering or changing values", {
  d <- qc_data()
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
})

test_that("report output has the stable compact contract", {
  r <- mad_qc(qc_data(), c(lower = "lower", upper = "upper"), output = "report", verbose = FALSE)
  expect_s3_class(r, "verymad_qc")
  expect_identical(names(r), c("flags", "thresholds", "settings"))
  expect_identical(names(r$flags), c("id", "lower", "upper", "mad_qc_outlier"))
  expect_identical(names(r$thresholds), c("metric", "direction", "transform", "median", "mad", "lower", "upper", "lower_raw", "upper_raw", "status"))
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

test_that("empty, missing, and insufficient metrics retain explicit status", {
  empty <- mad_qc(data.frame(a = numeric(), b = numeric()), c(a = "lower", b = "upper"), output = "report", verbose = FALSE)
  expect_equal(nrow(empty$flags), 0)
  expect_true(all(empty$thresholds$status == "all_missing"))
  expect_true(all(is.na(empty$thresholds[, c("median", "mad", "lower", "upper")])))

  missing <- mad_qc(data.frame(a = rep(NA_real_, 5)), c(a = "both"), output = "report", verbose = FALSE)
  expect_equal(missing$thresholds$status, "all_missing")
  expect_true(all(is.na(missing$flags$a)))
  expect_true(all(is.na(missing$flags$mad_qc_outlier)))

  expect_warning(insufficient <- mad_qc(data.frame(a = c(1, 2, NA, NA, NA), b = c(1, 2, 3, 4, 5)),
    c(a = "lower", b = "upper"), output = "report", min_n = 3, verbose = FALSE), "a")
  expect_equal(insufficient$thresholds$status, c("insufficient_n", "ok"))
})

test_that("observation identifiers fall back deterministically", {
  d <- data.frame(a = 1:5)
  attr(d, "row.names") <- c("dup", "dup", "", NA, "dup")
  r <- mad_qc(d, c(a = "lower"), output = "report", verbose = FALSE)
  expect_identical(r$flags$id, as.character(seq_len(nrow(d))))
})

test_that("zero MAD policies are explicit", {
  d <- data.frame(a = c(1, 2, 2, 2, 3))
  na <- mad_qc(d, c(a = "both"), output = "report", zero_mad = "na", verbose = FALSE)
  expect_equal(na$thresholds$status, "zero_mad"); expect_true(all(is.na(na$flags$a)))
  zero <- mad_qc(d, c(a = "both"), output = "report", zero_mad = "zero", verbose = FALSE)
  expect_equal(zero$thresholds$status, "zero_mad"); expect_true(all(is.logical(zero$flags$a)))
  expect_error(mad_qc(d, c(a = "both"), zero_mad = "error", verbose = FALSE), "a")
})

test_that("overwrite detects all conflicts before annotation", {
  d <- data.frame(a = 1:5, b = 1:5, a_mad_outlier = FALSE, b_mad_outlier = FALSE, mad_qc_outlier = FALSE)
  expect_error(mad_qc(d, c(a = "lower", b = "upper"), verbose = FALSE), "a_mad_outlier.*b_mad_outlier.*mad_qc_outlier")
  out <- mad_qc(d, c(a = "lower", b = "upper"), overwrite = TRUE, verbose = FALSE)
  expect_true(all(c("a_mad_outlier", "b_mad_outlier", "mad_qc_outlier") %in% names(out)))
})

test_that("mad_score matches manual median and MAD scaling", {
  x <- c(a = 10, b = 11, c = 11, d = 12, e = 12, f = 13, extreme = 100)
  expected <- (x - stats::median(x)) / stats::mad(x, constant = 1.4826)
  expect_equal(mad_score(x), expected)
  expect_type(mad_score(x), "double")
  expect_named(mad_score(x), names(x))
  expect_length(mad_score(x), length(x))
  expect_gt(abs(mad_score(x)[["extreme"]]), abs(as.numeric(scale(x))[[7]]))
})

test_that("mad_score handles missing values and zero MAD explicitly", {
  x <- c(a = 2, b = 2, missing = NA)
  expect_equal(mad_score(x), c(a = 0, b = 0, missing = NA))
  expect_true(all(is.na(mad_score(x, zero_mad = "na"))))
  expect_error(mad_score(x, zero_mad = "error"), "MAD is zero")
  expect_true(all(is.na(mad_score(c(1, NA, 3), na_rm = FALSE))))
})

test_that("limits and flags are direction-aware without inactive infinities", {
  x <- c(a = 0, b = 9, c = 10, d = 11, e = 20, missing = NA)
  lower <- mad_limits(x, nmads = 1, direction = "lower")
  upper <- mad_limits(x, nmads = 1, direction = "upper")
  both <- mad_limits(x, nmads = 1, direction = "both")
  expect_s3_class(lower, "tbl_df")
  expect_true(is.finite(lower$lower)); expect_true(is.na(lower$upper))
  expect_true(is.na(upper$lower)); expect_true(is.finite(upper$upper))
  expect_true(all(is.finite(c(both$lower, both$upper))))
  expect_false(any(is.infinite(unlist(list(lower, upper, both)))))
  expect_true(is_mad_outlier(x, 0.5, "lower")[["a"]])
  expect_false(is_mad_outlier(x, 0.5, "lower")[["e"]])
  expect_true(is_mad_outlier(x, 0.5, "upper")[["e"]])
  expect_true(is.na(is_mad_outlier(x)[["missing"]]))
  expect_equal(is_mad_outlier(x, 1, "both"), abs(mad_score(x)) > 1)
})

test_that("zero MAD flags follow all modes", {
  x <- c(2, 2, NA)
  expect_equal(is_mad_outlier(x), c(FALSE, FALSE, NA))
  expect_true(all(is.na(is_mad_outlier(x, zero_mad = "na"))))
  expect_error(is_mad_outlier(x, zero_mad = "error"), "MAD is zero")
})

test_that("vector validation rejects unsupported and non-finite inputs", {
  invalid <- list("x", factor("x"), matrix(1:4, 2), data.frame(x = 1:3),
    c(1, Inf), c(1, -Inf), c(1, NaN))
  for (x in invalid) expect_error(mad_score(x))
  expect_error(mad_score(1:3, constant = 0), "positive")
  expect_error(mad_limits(1:3, nmads = 0), "positive")
  expect_error(mad_limits(1:3, constant = Inf), "positive")
  expect_error(mad_limits(1:3, direction = "sideways"))
})

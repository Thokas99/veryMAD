test_that("is_outlier flags points beyond the threshold", {
  x <- c(1, 2, 2, 3, 100)
  expect_equal(is_outlier(x), c(FALSE, FALSE, FALSE, FALSE, TRUE))
})

test_that("is_outlier respects a custom threshold", {
  x <- c(1, 2, 2, 3, 100)
  expect_equal(
    is_outlier(x, threshold = 0.1),
    c(TRUE, FALSE, FALSE, TRUE, TRUE)
  )
})

test_that("is_outlier propagates missing values", {
  x <- c(1, 2, 2, 3, NA)
  expect_equal(is_outlier(x, na_rm = TRUE), c(FALSE, FALSE, FALSE, FALSE, NA))
})

test_that("is_outlier returns all NA when the scale is zero", {
  expect_equal(is_outlier(c(5, 5, 5, 5)), rep(NA, 4))
})

test_that("is_outlier errors on invalid threshold", {
  expect_snapshot(is_outlier(1:5, threshold = -1), error = TRUE)
  expect_snapshot(is_outlier(1:5, threshold = c(1, 2)), error = TRUE)
})

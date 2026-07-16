test_that("winsorize caps values beyond the threshold", {
  x <- c(1, 2, 2, 3, 100)
  centre <- median(x)
  scale <- mad(x, center = centre)
  upper <- centre + 3.5 * scale

  result <- winsorize(x)
  expect_equal(result[5], upper)
  expect_equal(result[1:4], x[1:4])
})

test_that("winsorize leaves x unchanged when the scale is zero", {
  x <- c(5, 5, 5, 5)
  expect_equal(winsorize(x), x)
})

test_that("winsorize handles missing values with na_rm", {
  x <- c(1, 2, 2, 3, NA)
  result <- winsorize(x, na_rm = TRUE)
  expect_equal(result[1:4], winsorize(c(1, 2, 2, 3)))
  expect_true(is.na(result[5]))
})

test_that("winsorize errors on invalid threshold", {
  expect_snapshot(winsorize(1:5, threshold = -1), error = TRUE)
  expect_snapshot(winsorize(1:5, threshold = c(1, 2)), error = TRUE)
})

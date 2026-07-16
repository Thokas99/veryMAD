test_that("mad_score returns 0 for the median value", {
  x <- c(1, 2, 2, 3, 100)
  expect_equal(mad_score(x)[2], 0)
})

test_that("mad_score matches manual calculation", {
  x <- c(1, 2, 2, 3, 100)
  centre <- median(x)
  scale <- mad(x, center = centre)
  expect_equal(mad_score(x), (x - centre) / scale)
})

test_that("mad_score handles missing values with na_rm", {
  x <- c(1, 2, 2, 3, NA)
  expect_equal(
    mad_score(x, na_rm = TRUE),
    c(mad_score(c(1, 2, 2, 3)), NA_real_)
  )
})

test_that("mad_score returns all NA when the scale is zero", {
  expect_equal(mad_score(c(5, 5, 5, 5)), rep(NA_real_, 4))
})

test_that("mad_score errors on non-numeric input", {
  expect_snapshot(mad_score("a"), error = TRUE)
})

test_that("mad_score errors on invalid constant", {
  expect_snapshot(mad_score(1:5, constant = -1), error = TRUE)
  expect_snapshot(mad_score(1:5, constant = c(1, 2)), error = TRUE)
})

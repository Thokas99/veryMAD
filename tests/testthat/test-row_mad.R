test_that("row_mad matches stats::mad() per row", {
  m <- matrix(c(1, 2, 2, 3, 100, 10, 20, 20, 30, 40), nrow = 2, byrow = TRUE)
  result <- row_mad(m)

  expect_equal(unname(result), c(mad(m[1, ]), mad(m[2, ])))
})

test_that("row_mad can compute across columns", {
  m <- matrix(c(1, 2, 2, 3, 100, 10, 20, 20, 30, 40), nrow = 2, byrow = TRUE)
  result <- row_mad(m, margin = 2)

  expect_equal(unname(result), apply(m, 2, mad))
})

test_that("row_mad preserves names", {
  m <- matrix(1:4, nrow = 2, dimnames = list(c("g1", "g2"), c("s1", "s2")))
  expect_equal(names(row_mad(m)), c("g1", "g2"))
  expect_equal(names(row_mad(m, margin = 2)), c("s1", "s2"))
})

test_that("row_mad errors on non-matrix input", {
  expect_snapshot(row_mad(data.frame(x = 1:5)), error = TRUE)
})

test_that("row_mad errors on invalid margin", {
  m <- matrix(1:4, nrow = 2)
  expect_snapshot(row_mad(m, margin = 3), error = TRUE)
})

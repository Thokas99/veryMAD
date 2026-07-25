test_that("matrix helpers preserve names and orientation", {
  x <- matrix(c(1, 2, 3, 2, 2, 8), nrow = 2, byrow = TRUE,
              dimnames = list(c("g1", "g2"), c("s1", "s2", "s3")))
  expect_named(row_mad(x), rownames(x))
  expect_named(col_mad(x), colnames(x))
  rows <- robust_scale(x, "rows")
  columns <- robust_scale(x, "columns")
  expect_equal(dim(rows), dim(x)); expect_equal(dimnames(rows), dimnames(x))
  expect_equal(dim(columns), dim(x)); expect_equal(dimnames(columns), dimnames(x))
})

test_that("matrix missing and zero MAD behavior is explicit", {
  x <- matrix(c(1, 1, NA, 1, 2, 3), nrow = 2, byrow = TRUE)
  expect_equal(robust_scale(x)[1, ], c(0, 0, NA))
  expect_true(all(is.na(robust_scale(x, zero_mad = "na")[1, ])))
  expect_error(robust_scale(x, zero_mad = "error"), "MAD is zero")
  expect_error(robust_scale(x, "cells"))
  expect_error(row_mad(data.frame(x = 1:3)))
  expect_error(row_mad(matrix(c(1, Inf), nrow = 1)))
  sparse <- structure(matrix(1, 1, 1), class = c("dgCMatrix", "sparseMatrix"))
  expect_error(row_mad(sparse), "Sparse matrices")
})

test_that("top MAD features ranks rows only", {
  x <- rbind(stable = c(1, 1, 1), variable = c(1, 2, 10))
  expect_equal(top_mad_features(x, 1), "variable")
})

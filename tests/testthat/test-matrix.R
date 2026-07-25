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
  sparse <- structure(matrix(1, 1, 1), class = "sparseMatrix")
  expect_error(row_mad(sparse), "dgCMatrix")
})

test_that("top MAD features ranks rows only", {
  x <- rbind(stable = c(1, 1, 1), variable = c(1, 2, 10))
  expect_equal(top_mad_features(x, 1), "variable")
})

test_that("dgCMatrix MADs match dense calculations without densifying", {
  skip_if_not_installed("Matrix")
  skip_if_not_installed("sparseMatrixStats")
  dense <- rbind(zero = c(0, 0, 0, 0), tied_a = c(0, 1, 0, 1),
                 tied_b = c(0, 1, 0, 1), missing = c(0, NA, 2, 0))
  colnames(dense) <- paste0("s", 1:4)
  sparse <- Matrix::Matrix(dense, sparse = TRUE)
  expect_s4_class(sparse, "dgCMatrix")
  expect_equal(row_mad(sparse, na_rm = TRUE), row_mad(dense, na_rm = TRUE))
  expect_equal(col_mad(sparse, na_rm = TRUE), col_mad(dense, na_rm = TRUE))
  expect_named(row_mad(sparse), rownames(dense))
  expect_named(col_mad(sparse), colnames(dense))
  expect_equal(top_mad_features(sparse, 4), top_mad_features(dense, 4))
  expect_equal(top_mad_features(sparse, 4)[1:2], c("tied_a", "tied_b"))
  expect_error(robust_scale(sparse), "destroy sparsity")
})

test_that("sparse MAD dependency failure is informative", {
  skip_if_not_installed("Matrix")
  sparse <- Matrix::sparseMatrix(i = 1:2, j = 1:2, x = c(1, 1))
  testthat::local_mocked_bindings(
    .require_namespace = function(package, feature) stop("sparseMatrixStats unavailable", call. = FALSE),
    .package = "veryMAD"
  )
  expect_error(row_mad(sparse), "sparseMatrixStats unavailable")
})

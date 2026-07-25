test_that("mad_scale vectors equal mad_score", {
  x <- c(a = 1, b = 2, c = 10, missing = NA)
  expect_equal(mad_scale(x), mad_score(x))
  expect_named(mad_scale(x), names(x))
})

test_that("mad_scale scales rows and columns like separate vectors", {
  x <- matrix(c(1, 2, 7, 2, 4, 10, 3, 8, 20), nrow = 3, byrow = TRUE,
    dimnames = list(c("g1", "g2", "g3"), c("s1", "s2", "s3")))
  rows <- mad_scale(x, "rows")
  columns <- mad_scale(x, "columns")
  expected_rows <- do.call(rbind, lapply(seq_len(nrow(x)), function(i) mad_score(x[i, ])))
  expected_columns <- do.call(cbind, lapply(seq_len(ncol(x)), function(j) mad_score(x[, j])))
  dimnames(expected_rows) <- dimnames(x); dimnames(expected_columns) <- dimnames(x)
  expect_equal(rows, expected_rows)
  expect_equal(columns, expected_columns)
  expect_equal(dim(rows), dim(x)); expect_equal(dimnames(rows), dimnames(x))
  expect_equal(dim(columns), dim(x)); expect_equal(dimnames(columns), dimnames(x))
  expect_true(is.matrix(rows)); expect_type(rows, "double")
})

test_that("mad_scale supports integer, one-row, one-column, and missing inputs", {
  integer_matrix <- matrix(1:6, nrow = 2)
  expect_type(mad_scale(integer_matrix), "double")
  expect_equal(dim(mad_scale(matrix(c(1, 2, 4), nrow = 1))), c(1L, 3L))
  expect_equal(dim(mad_scale(matrix(c(1, 2, 4), ncol = 1), "columns")), c(3L, 1L))
  x <- matrix(c(1, NA, 3, 2, 4, 8), nrow = 2, byrow = TRUE)
  expect_true(is.na(mad_scale(x)[1, 2]))
  expect_true(all(is.na(mad_scale(x, na_rm = FALSE)[1, ])))
})

test_that("mad_scale handles zero and missing spreads", {
  x <- rbind(TP53 = c(2, 2, NA), variable = c(1, 2, 4), missing = c(NA, NA, NA))
  zero <- mad_scale(x, zero_mad = "zero")
  expect_equal(zero["TP53", ], c(0, 0, NA))
  expect_true(all(is.na(zero["missing", ])))
  expect_true(all(is.na(mad_scale(x, zero_mad = "na")["TP53", ])))
  expect_error(mad_scale(x, zero_mad = "error"), "row `TP53`")
  y <- cbind(Sample_01 = c(1, 2, 3), Sample_03 = c(5, 5, 5))
  expect_error(mad_scale(y, "columns", zero_mad = "error"), "column `Sample_03`")
})

test_that("mad_scale rejects unsupported inputs and sparse matrices", {
  expect_error(mad_scale(data.frame(x = 1:3)), "ordinary numeric matrix")
  expect_error(mad_scale(tibble::tibble(x = 1:3)), "ordinary numeric matrix")
  expect_error(mad_scale(array(1:8, c(2, 2, 2))), "ordinary numeric matrix")
  expect_error(mad_scale(matrix(letters[1:4], 2)), "numeric matrix")
  expect_error(mad_scale(matrix(c(1, Inf), 1)), "Inf")
  sparse <- structure(matrix(1, 1, 1), class = c("dgCMatrix", "sparseMatrix"))
  expect_error(mad_scale(sparse), "ordinary dense numeric matrices only")
  expect_error(mad_scale(matrix(1:4, 2), margin = "cells"))
})

test_that("mad_scale uses matrixStats margin functions", {
  code <- paste(deparse(body(mad_scale)), collapse = "\n")
  expect_match(code, "matrixStats::rowMedians", fixed = TRUE)
  expect_match(code, "matrixStats::rowMads", fixed = TRUE)
  expect_match(code, "matrixStats::colMedians", fixed = TRUE)
  expect_match(code, "matrixStats::colMads", fixed = TRUE)
  expect_false(grepl("apply\\(", code))
})

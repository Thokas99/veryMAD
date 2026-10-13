test_that("mad_scale handles numeric vectors and names", {
  x <- c(a = 1, b = 2, c = 100)
  out <- mad_scale(x)
  expect_type(out, "double"); expect_identical(names(out), names(x))
  expect_equal(mad_scale(x, center = FALSE, scale = FALSE), x)
  expect_equal(mad_scale(numeric()), numeric())
  expect_error(mad_scale(c(1, 2, Inf)), "Inf")
  expect_error(mad_scale(c(1, 2, -Inf)), "Inf")
  expect_error(mad_scale(c(1, 2, NaN)), "NaN")
  expect_error(mad_scale(letters), "numeric")
})

test_that("mad_scale handles matrices, margins, and data frames", {
  m <- matrix(1:12, nrow = 3, dimnames = list(paste0("r", 1:3), paste0("c", 1:4)))
  rows <- mad_scale(m, margin = 1)
  columns <- mad_scale(m, margin = 2)
  expect_true(is.matrix(rows)); expect_identical(dim(rows), dim(m)); expect_identical(dimnames(rows), dimnames(m))
  expect_equal(rows, mad_scale(m, margin = "rows")); expect_equal(columns, mad_scale(m, margin = "columns"))
  expect_identical(dim(mad_scale(matrix(numeric(), nrow = 0, ncol = 2))), c(0L, 2L))
  expect_true(is.matrix(mad_scale(as.data.frame(m))))
  expect_error(mad_scale(data.frame(a = 1:3, b = letters[1:3])), "numeric")
  expect_error(mad_scale(m, margin = 3), "margin")
})

test_that("mad_scale handles missingness and all zero-MAD policies", {
  x <- c(1, 2, 3, NA_real_)
  expect_true(is.na(mad_scale(x, na_rm = FALSE)[4]))
  expect_true(is.na(mad_scale(matrix(c(1, 2, NA, 4), nrow = 2), margin = 1, na_rm = FALSE)[1, 1]))
  expect_true(all(is.na(mad_scale(c(1, 1, 1), zero_mad = "na"))))
  expect_equal(mad_scale(c(1, 1, 1), zero_mad = "zero"), c(0, 0, 0))
  expect_error(mad_scale(c(1, 1, 1), zero_mad = "error"), "zero")
  expect_true(all(is.na(mad_scale(matrix(1, nrow = 2, ncol = 2), zero_mad = "na"))))
})

test_that("mad_scale defaults to unresolved zero MAD and preserves vector names", {
  expect_identical(eval(formals(mad_scale)$zero_mad), c("na", "zero", "error"))

  constant <- c(a = 4, b = 4, c = 4)
  expect_true(all(is.na(mad_scale(constant))))
  expect_identical(names(mad_scale(constant)), names(constant))

  mostly_constant <- c(a = 1, b = 1, c = 1, outlier = 100)
  expect_true(all(is.na(mad_scale(mostly_constant))))
  expect_equal(mad_scale(mostly_constant, zero_mad = "zero"),
               c(a = 0, b = 0, c = 0, outlier = 99))
  expect_equal(mad_scale(mostly_constant, center = FALSE, zero_mad = "zero"),
               mostly_constant)
})

test_that("zero-MAD handling is skipped when scaling is disabled", {
  x <- c(a = 1, b = 1, c = 1, outlier = 100)
  expect_equal(mad_scale(x, scale = FALSE, zero_mad = "error"),
               c(a = 0, b = 0, c = 0, outlier = 99))
  expect_equal(mad_scale(x, center = FALSE, scale = FALSE, zero_mad = "error"), x)

  m <- rbind(constant = c(1, 1, 1), variable = c(1, 2, 4))
  expect_silent(mad_scale(m, margin = 1, scale = FALSE, zero_mad = "error"))
  expect_silent(mad_scale(m, margin = 1, center = FALSE, scale = FALSE,
                          zero_mad = "error"))
})

test_that("zero-MAD matrix rows and columns preserve dimensions and dimnames", {
  rows <- rbind(
    constant = c(2, 2, 2, 2),
    variable = c(1, 2, 4, 8)
  )
  dimnames(rows) <- list(c("constant", "variable"), paste0("cell", 1:4))
  row_scaled <- mad_scale(rows, margin = 1)
  expect_identical(dim(row_scaled), dim(rows))
  expect_identical(dimnames(row_scaled), dimnames(rows))
  expect_true(all(is.na(row_scaled["constant", ])))
  expect_equal(unname(mad_scale(rows, margin = 1, zero_mad = "zero")["constant", ]),
               rep(0, ncol(rows)))
  expect_equal(mad_scale(rows, margin = 1, center = FALSE,
                         zero_mad = "zero")["constant", ],
               rows["constant", ])

  columns <- cbind(
    constant = c(3, 3, 3, 3),
    variable = c(1, 2, 4, 8)
  )
  rownames(columns) <- paste0("cell", 1:4)
  column_scaled <- mad_scale(columns, margin = 2)
  expect_identical(dim(column_scaled), dim(columns))
  expect_identical(dimnames(column_scaled), dimnames(columns))
  expect_true(all(is.na(column_scaled[, "constant"])))
  expect_equal(unname(mad_scale(columns, margin = 2, zero_mad = "zero")[, "constant"]),
               rep(0, nrow(columns)))

  expect_error(mad_scale(rows, margin = 1, zero_mad = "error"), "row 1")
  expect_error(mad_scale(columns, margin = 2, zero_mad = "error"), "column 1")
})

test_that("missing vectors and matrix margins retain explicit missingness", {
  all_missing <- c(a = NA_real_, b = NA_real_)
  expect_true(all(is.na(mad_scale(all_missing))))
  expect_identical(names(mad_scale(all_missing)), names(all_missing))

  mixed <- c(a = 1, b = NA_real_, c = 2, d = 3)
  scaled <- mad_scale(mixed)
  expect_true(is.na(scaled["b"]))
  expect_true(all(is.finite(scaled[c("a", "c", "d")])) )

  m <- rbind(constant_with_missing = c(1, 1, NA_real_, 1),
             variable_with_missing = c(1, 2, NA_real_, 8))
  out <- mad_scale(m, margin = 1)
  expect_true(all(is.na(out["constant_with_missing", ])))
  expect_true(is.na(out["variable_with_missing", 3]))
  expect_true(all(is.finite(out["variable_with_missing", -3])))
})

test_that("mad_scale agrees with stats MAD", {
  x <- c(a = 1, b = 2, c = 4, d = 5, e = 9, f = 20, missing = NA_real_)
  centre <- stats::median(x, na.rm = TRUE)
  spread <- stats::mad(x, center = centre, constant = 1.4826, na.rm = TRUE)
  expect_equal(mad_scale(x), (x - centre) / spread)
  expect_equal(mad_scale(x, center = FALSE, scale = TRUE), x / spread)
})

test_that("mad_qc matches an equivalent scuttle outlier calculation", {
  skip_if_not_installed("scuttle")
  depth <- c(2, 3, 4, 5, 6, 8, 10, 30)
  report <- mad_qc(
    data.frame(depth = depth), c(depth = "upper"), nmads = 2,
    transform = c(depth = "log1p"), min_n = 5, output = "report",
    verbose = FALSE
  )
  reference <- scuttle::isOutlier(
    log1p(depth), nmads = 2, type = "higher", log = FALSE,
    min.diff = NA_real_
  )
  expect_identical(report$flags$depth, as.logical(reference))
  transformed <- log1p(depth[is.finite(depth)])
  expect_equal(report$thresholds$median, stats::median(transformed))
  expect_equal(report$thresholds$mad,
               stats::mad(transformed, center = stats::median(transformed), constant = 1.4826))

  constant <- rep(4, 6)
  zero_mad_report <- mad_qc(
    data.frame(value = constant), c(value = "both"), zero_mad = "zero",
    output = "report", verbose = FALSE
  )
  zero_mad_reference <- scuttle::isOutlier(
    constant, nmads = 3, type = "both", log = FALSE, min.diff = NA_real_
  )
  expect_identical(zero_mad_report$flags$value, as.logical(zero_mad_reference))
})

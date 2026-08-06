test_that("mad_qc annotates explicit metrics without filtering", {
  d <- data.frame(a = c(1, 10, 11, 12, 100), b = c(1, 2, 3, 4, 20), row.names = letters[1:5])
  out <- mad_qc(d, c(a = "lower", b = "upper"), verbose = FALSE)
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), nrow(d)); expect_identical(rownames(out), rownames(d))
  expect_true(all(c("a_mad_outlier", "b_mad_outlier", "mad_qc_outlier") %in% names(out)))
  expect_identical(d$a, out$a); expect_true(any(out$mad_qc_outlier))
})

test_that("mad_qc reports compact thresholds and three-valued flags", {
  d <- data.frame(a = c(NA, 1, 2, 3, 4), b = c(1, 1, 1, 1, 1), row.names = letters[1:5])
  expect_warning(r <- mad_qc(d, c(a = "both", b = "upper"), output = "report", verbose = FALSE), "Insufficient")
  expect_s3_class(r, "verymad_qc"); expect_named(r, c("flags", "thresholds", "settings"))
  expect_equal(names(r$flags), c("id", "a", "b", "mad_qc_outlier"))
  expect_equal(names(r$thresholds), c("metric", "direction", "transform", "median", "mad", "lower", "upper", "lower_raw", "upper_raw", "status"))
  expect_true(all(r$thresholds$status %in% c("ok", "insufficient_n", "all_missing", "zero_mad")))
  expect_true(all(is.na(r$flags$b))); expect_true(is.na(r$flags$mad_qc_outlier[1]))
})

test_that("mad_qc validates metrics, transforms, and domains", {
  d <- data.frame(a = 1:5, b = letters[1:5])
  expect_silent(mad_qc(d, c(a = "lower"), verbose = FALSE))
  expect_error(mad_qc(d, c(a = "lower", a = "upper"), verbose = FALSE), "unique")
  expect_error(mad_qc(d, c(b = "lower"), verbose = FALSE), "numeric")
  expect_error(mad_qc(d, c(nope = "lower"), verbose = FALSE), "Missing")
  expect_error(mad_qc(d, c(a = "sideways"), verbose = FALSE), "directions")
  expect_error(mad_qc(data.frame(a = -1:3), c(a = "lower"), transform = "log1p", verbose = FALSE), "non-negative")
  expect_error(mad_qc(data.frame(a = 0:4), c(a = "lower"), transform = "log10", verbose = FALSE), "positive")
})

test_that("mad_qc handles zero MAD policies and partial transforms", {
  d <- data.frame(a = c(1, 2, 2, 2, 3), b = 1:5)
  r <- mad_qc(d, c(a = "both", b = "lower"), transform = c(a = "log1p"), output = "report", zero_mad = "na", verbose = FALSE)
  expect_equal(r$settings$transform, c(a = "log1p", b = "none"))
  z <- mad_qc(data.frame(a = rep(2, 5)), c(a = "both"), output = "report", zero_mad = "zero", verbose = FALSE)
  expect_equal(z$thresholds$status, "zero_mad"); expect_false(any(z$flags$a, na.rm = TRUE))
  expect_error(mad_qc(data.frame(a = rep(2, 5)), c(a = "both"), zero_mad = "error", verbose = FALSE), "a")
})

test_that("matrix and tibble-shaped inputs work", {
  x <- matrix(1:15, nrow = 5, dimnames = list(letters[1:5], c("a", "b", "c")))
  out <- mad_qc(x, c(a = "lower", b = "upper"), verbose = FALSE)
  expect_true(is.data.frame(out)); expect_equal(nrow(out), 5)
  tibble <- as.data.frame(x); class(tibble) <- c("tbl_df", "tbl", "data.frame")
  expect_s3_class(mad_qc(tibble, c(a = "lower"), verbose = FALSE), "data.frame")
})

test_that("Seurat input uses the same engine when available", {
  skip_if_not_installed("SeuratObject")
  counts <- matrix(1:25, nrow = 5, dimnames = list(paste0("g", 1:5), paste0("c", 1:5)))
  object <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  object[["qc"]] <- c(1, 2, 3, 4, 100)
  out <- mad_qc(object, c(qc = "upper"), verbose = FALSE)
  expect_true(inherits(out, "Seurat"))
  expect_true(all(c("qc_mad_outlier", "mad_qc_outlier") %in% names(out[[]])))
})

test_that("mad_scale supports vectors, matrices, names, missingness, and alias", {
  x <- c(a = 1, b = 2, c = 100)
  expect_equal(mad_z_score(x), mad_scale(x))
  expect_equal(names(mad_scale(x)), names(x))
  m <- matrix(1:12, nrow = 3, dimnames = list(paste0("r", 1:3), paste0("c", 1:4)))
  expect_equal(dim(mad_scale(m, margin = 1)), dim(m)); expect_equal(dimnames(mad_scale(m)), dimnames(m))
  expect_true(is.matrix(mad_scale(as.data.frame(m))))
  expect_true(all(is.na(mad_scale(matrix(1, nrow = 2, ncol = 2), zero_mad = "na"))))
  expect_error(mad_scale(c(1, 1, 1), zero_mad = "error"), "zero")
})

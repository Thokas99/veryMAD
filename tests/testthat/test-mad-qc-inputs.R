test_that("mad_qc validates the explicit metric contract", {
  d <- data.frame(a = 1:5, b = letters[1:5])
  expect_silent(mad_qc(d, c(a = "lower"), verbose = FALSE))
  expect_error(mad_qc(d, c(a = "lower", a = "upper"), verbose = FALSE), "unique")
  expect_error(mad_qc(d, c(b = "lower"), verbose = FALSE), "numeric")
  expect_error(mad_qc(d, c(nope = "lower"), verbose = FALSE), "Missing")
  expect_error(mad_qc(d, c(a = "sideways"), verbose = FALSE), "directions")
  expect_error(mad_qc(d, c(a = "lower"), nmads = 0, verbose = FALSE), "positive")
  expect_error(mad_qc(d, c(a = "lower"), min_n = 0, verbose = FALSE), "positive")
  expect_error(mad_qc(d, c(a = "lower"), min_n = 1.5, verbose = FALSE), "integer")
  expect_error(mad_qc(d, c(a = "lower"), transform = "unknown", verbose = FALSE), "Transformations")
  expect_error(mad_qc(d, c(a = "lower"), transform = c(unknown = "log1p"), verbose = FALSE), "selected metrics")
  expect_false("combine" %in% names(formals(mad_qc)))
})

test_that("reserved report and annotation names cannot be selected as metrics", {
  for (reserved in c("id", "mad_qc_outlier")) {
    d <- data.frame(value = 1:5)
    names(d) <- reserved
    metrics <- stats::setNames("upper", reserved)

    expect_error(mad_qc(d, metrics, verbose = FALSE), "[Rr]eserved.*(id|mad_qc_outlier)")
    expect_error(mad_qc(d, metrics, output = "report", verbose = FALSE),
                 "[Rr]eserved.*(id|mad_qc_outlier)")
  }
})

test_that("mad_qc validates transformations and finite values", {
  expect_error(mad_qc(data.frame(a = -1:3), c(a = "lower"), transform = "log1p", verbose = FALSE), "non-negative")
  expect_error(mad_qc(data.frame(a = 0:4), c(a = "lower"), transform = "log10", verbose = FALSE), "positive")
  expect_error(mad_qc(data.frame(a = c(1:4, Inf)), c(a = "lower"), verbose = FALSE), "non-finite")
  expect_error(mad_qc(data.frame(a = c(1:4, NaN)), c(a = "lower"), verbose = FALSE), "non-finite")
  expect_error(mad_qc(data.frame(a = c(1:4, -Inf)), c(a = "lower"), verbose = FALSE), "non-finite")
})

test_that("matrix, tibble-like, and optional Seurat inputs work", {
  x <- matrix(1:15, nrow = 5, dimnames = list(letters[1:5], c("a", "b", "c")))
  out <- mad_qc(x, c(a = "lower", b = "upper"), verbose = FALSE)
  expect_true(is.data.frame(out)); expect_equal(nrow(out), 5)
  report <- mad_qc(x, c(a = "lower"), output = "report", verbose = FALSE)
  expect_identical(report$flags$id, rownames(x))
  expect_equal(nrow(report$flags), nrow(x))
  tibble <- as.data.frame(x); class(tibble) <- c("tbl_df", "tbl", "data.frame")
  expect_s3_class(mad_qc(tibble, c(a = "lower"), verbose = FALSE), "data.frame")

  skip_if_not_installed("SeuratObject")
  assay_counts <- function(x) {
    if (utils::packageVersion("SeuratObject") >= "5.0.0") {
      SeuratObject::LayerData(x, assay = "RNA", layer = "counts")
    } else {
      SeuratObject::GetAssayData(x, assay = "RNA", slot = "counts")
    }
  }
  counts <- matrix(1:25, nrow = 5, dimnames = list(paste0("g", 1:5), paste0("c", 1:5)))
  object <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  object[["qc"]] <- c(1, 2, 3, 4, 100)
  counts_before <- assay_counts(object)
  report <- mad_qc(object, c(qc = "upper"), output = "report", verbose = FALSE)
  out <- mad_qc(object, c(qc = "upper"), verbose = FALSE)
  expect_true(inherits(out, "Seurat"))
  expect_true(all(c("qc_mad_outlier", "mad_qc_outlier") %in% names(out[[]])))
  expect_identical(report$flags$id, colnames(object))
  expect_identical(rownames(out[[]]), colnames(out))
  expect_identical(out[[]]$qc_mad_outlier, report$flags$qc)
  expect_identical(assay_counts(out), counts_before)
})

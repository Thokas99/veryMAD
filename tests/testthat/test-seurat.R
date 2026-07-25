test_that("Seurat integration reports and annotates by cell name", {
  skip_if_not_installed("SeuratObject")
  counts <- matrix(c(1, 0, 2, 0, 3, 10), nrow = 2,
                   dimnames = list(c("g1", "g2"), c("c1", "c2", "c3")))
  object <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  object[["percent.mt"]] <- c(1, 2, 50)
  original <- object
  report <- mad_qc_seurat(object, metrics = c(nCount_RNA = "upper", percent.mt = "upper"), action = "report")
  expect_s3_class(report, "mad_qc")
  expect_equal(report$id[report$metric == "percent.mt"], colnames(object))
  annotated <- mad_qc_seurat(object, metrics = c(percent.mt = "upper"), nmads = 1)
  expect_true(all(c("percent.mt_mad_outlier", "mad_qc_outlier") %in% colnames(annotated[[]])))
  expect_false("mad_qc_outlier" %in% colnames(original[[]]))
  expect_error(mad_qc_seurat(object, metrics = c(nope = "upper"), action = "report"), "Missing metric")
})

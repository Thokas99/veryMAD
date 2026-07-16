skip_if_not_installed("SeuratObject")

test_that("mad_qc dispatches on a Seurat object and writes @meta.data", {
  counts <- matrix(
    1:20,
    nrow = 4,
    dimnames = list(paste0("gene", 1:4), paste0("cell", 1:5))
  )
  obj <- SeuratObject::CreateSeuratObject(counts = counts)
  obj$percent.mt <- c(2, 3, 2.5, 3, 40)

  result <- mad_qc(
    obj,
    metrics = c(nCount_RNA = "lower", percent.mt = "upper")
  )

  expect_s4_class(result, "Seurat")
  expect_equal(
    result@meta.data$percent.mt_outlier,
    c(FALSE, FALSE, FALSE, FALSE, TRUE)
  )
  expect_true("nCount_RNA_outlier" %in% names(result@meta.data))
  expect_equal(
    result@meta.data$mad_qc_outlier,
    result@meta.data$nCount_RNA_outlier | result@meta.data$percent.mt_outlier
  )
})

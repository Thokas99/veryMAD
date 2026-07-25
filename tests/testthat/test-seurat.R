seurat_fixture <- function() {
  counts <- matrix(c(1, 0, 2, 0, 3, 10), nrow = 2,
    dimnames = list(c("g1", "g2"), c("c1", "c2", "c3")))
  object <- suppressWarnings(SeuratObject::CreateSeuratObject(counts))
  object[["percent.mt"]] <- c(1, 50, NA)
  object
}

test_that("Seurat report mode is default and leaves the object unchanged", {
  skip_if_not_installed("SeuratObject")
  object <- seurat_fixture(); original <- object
  report <- mad_qc_seurat(object, metrics = c(percent.mt = "upper"), nmads = 1)
  expect_s3_class(report, "mad_qc")
  expect_s3_class(report, "tbl_df")
  expect_null(attr(report, "group_by"))
  expect_equal(report$id, colnames(object))
  expect_identical(object, original)
})

test_that("Seurat annotation adds only aligned three-valued flags", {
  skip_if_not_installed("SeuratObject")
  object <- seurat_fixture(); original_metadata <- object[[]]
  original_cells <- colnames(object)
  original_assays <- SeuratObject::Assays(object)
  original_active <- SeuratObject::DefaultAssay(object)
  original_idents <- SeuratObject::Idents(object)
  original_reductions <- SeuratObject::Reductions(object)
  annotated <- mad_qc_seurat(object, metrics = c(percent.mt = "upper"), nmads = 0.5,
    action = "annotate")
  added <- setdiff(colnames(annotated[[]]), colnames(original_metadata))
  expect_equal(added, c("percent.mt_mad_outlier", "mad_qc_outlier"))
  expect_equal(colnames(annotated), original_cells)
  expect_equal(ncol(annotated), ncol(object))
  expect_equal(SeuratObject::Assays(annotated), original_assays)
  expect_equal(SeuratObject::DefaultAssay(annotated), original_active)
  expect_equal(SeuratObject::Idents(annotated), original_idents)
  expect_equal(SeuratObject::Reductions(annotated), original_reductions)
  expect_equal(annotated[[]][names(original_metadata)], original_metadata)
  expect_equal(annotated[[]]$mad_qc_outlier, c(FALSE, TRUE, NA))
  expect_false("mad_qc_outlier" %in% colnames(object[[]]))
})

test_that("Seurat annotation protects flag collisions", {
  skip_if_not_installed("SeuratObject")
  object <- seurat_fixture()
  annotated <- mad_qc_seurat(object, metrics = c(percent.mt = "upper"), action = "annotate")
  expect_error(mad_qc_seurat(annotated, metrics = c(percent.mt = "upper"), action = "annotate"), "already exist")
  overwritten <- mad_qc_seurat(annotated, metrics = c(percent.mt = "upper"),
    action = "annotate", overwrite = TRUE)
  expect_equal(colnames(overwritten), colnames(object))
  expect_error(mad_qc_seurat(object, metrics = c(nope = "upper")), "Missing metric")
})

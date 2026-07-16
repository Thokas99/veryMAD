test_that("mad_qc flags observations beyond the threshold", {
  cell_metadata <- data.frame(
    nCount_RNA = c(500, 600, 550, 20000, 580),
    percent.mt = c(2, 3, 2.5, 3, 40)
  )
  qc <- mad_qc(
    cell_metadata,
    metrics = c(nCount_RNA = "lower", percent.mt = "upper"),
    nmads = 3
  )

  expect_setequal(qc$metric, c("nCount_RNA", "percent.mt"))
  expect_equal(qc$is_outlier[qc$metric == "percent.mt"], c(rep(FALSE, 4), TRUE))
})

test_that("mad_qc applies log10 transforms", {
  cell_metadata <- data.frame(nCount_RNA = c(500, 600, 550, 20000, 580))
  qc <- mad_qc(
    cell_metadata,
    metrics = c(nCount_RNA = "lower"),
    transform = c(nCount_RNA = "log10")
  )

  expect_equal(qc$value, log10(cell_metadata$nCount_RNA))
})

test_that("mad_qc computes thresholds within groups", {
  cell_metadata <- data.frame(
    nCount_RNA = c(500, 600, 550, 5000, 6000, 5500),
    sample = rep(c("a", "b"), each = 3)
  )
  qc <- mad_qc(
    cell_metadata,
    metrics = c(nCount_RNA = "both"),
    group_by = "sample"
  )

  expect_equal(length(unique(qc$median)), 2)
  expect_false(any(qc$is_outlier))
})

test_that("mad_qc errors on invalid inputs", {
  cell_metadata <- data.frame(nCount_RNA = c(500, 600, 550))
  expect_snapshot(mad_qc(cell_metadata, metrics = c(nCount_RNA = "sideways")), error = TRUE)
  expect_snapshot(mad_qc(cell_metadata, metrics = c("nCount_RNA")), error = TRUE)
  expect_snapshot(mad_qc(cell_metadata, metrics = c(missing_col = "lower")), error = TRUE)
})

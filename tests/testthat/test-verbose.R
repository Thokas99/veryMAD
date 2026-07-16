test_that("verbose defaults stay quiet", {
  expect_silent(mad_score(1:5))
  expect_silent(is_outlier(1:5))
  expect_silent(winsorize(1:5))
  expect_silent(row_mad(matrix(1:9, nrow = 3)))
  expect_silent(robust_scale(matrix(1:9, nrow = 3)))
  expect_silent(select_variable_features(matrix(1:9, nrow = 3), n = 1))
  expect_silent(mad_qc(data.frame(x = 1:5), metrics = c(x = "both")))
})

test_that("verbose output reports useful summaries", {
  expect_message(mad_score(1:5, verbose = TRUE), "Calculated MAD scores")
  expect_message(is_outlier(c(1, 2, 100), verbose = TRUE), "Flagged")
  expect_message(winsorize(c(1, 2, 100), verbose = TRUE), "Winsorized")
  expect_message(row_mad(matrix(1:9, nrow = 3), verbose = TRUE), "Calculated MAD")
  expect_message(robust_scale(matrix(1:9, nrow = 3), verbose = TRUE), "Robust-scaled")
  expect_message(select_variable_features(matrix(1:9, nrow = 3), n = 1, verbose = TRUE), "Selected")
  expect_message(
    mad_qc(data.frame(x = 1:5), metrics = c(x = "both"), verbose = TRUE),
    "QC flagged"
  )
})

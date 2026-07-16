test_that("count and threshold arguments reject non-finite values", {
  x <- matrix(1:9, nrow = 3)

  expect_error(select_variable_features(x, n = 1.5), "finite positive whole")
  expect_error(select_variable_features(x, n = Inf), "finite positive whole")
  expect_error(is_outlier(1:3, threshold = Inf), "finite positive")
  expect_error(winsorize(1:3, threshold = Inf), "finite positive")
  expect_error(mad_qc(data.frame(x = 1:3), metrics = c(x = "both"), nmads = Inf), "finite positive")
})

test_that("select_variable_features picks the top n rows by MAD", {
  m <- matrix(
    c(1, 50, 100, 150, 200, 10, 20, 20, 30, 40, 1, 1, 1, 1, 1),
    nrow = 3,
    byrow = TRUE
  )
  rownames(m) <- c("high_var", "mid_var", "no_var")

  expect_equal(select_variable_features(m, n = 1), "high_var")
  expect_equal(
    select_variable_features(m, n = 2),
    c("high_var", "mid_var")
  )
})

test_that("select_variable_features can rank columns", {
  m <- matrix(c(1, 2, 2, 3, 100, 10, 20, 20, 30, 40), nrow = 2, byrow = TRUE)
  colnames(m) <- letters[1:5]

  expect_equal(select_variable_features(m, n = 1, margin = 2), "e")
})

test_that("select_variable_features returns indices when unnamed", {
  m <- matrix(c(1, 2, 2, 3, 100, 10, 20, 20, 30, 40), nrow = 2, byrow = TRUE)
  expect_equal(select_variable_features(m, n = 1), 2)
})

test_that("select_variable_features caps n at the number of features", {
  m <- matrix(c(1, 2, 2, 3, 100, 10, 20, 20, 30, 40), nrow = 2, byrow = TRUE)
  expect_length(select_variable_features(m, n = 100), 2)
})

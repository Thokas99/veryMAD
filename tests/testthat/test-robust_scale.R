test_that("robust_scale median-centers and MAD-scales rows by default", {
  m <- matrix(c(1, 2, 2, 3, 100, 10, 20, 20, 30, 40), nrow = 2, byrow = TRUE)
  result <- robust_scale(m)

  expect_equal(result[1, ], mad_score(m[1, ]))
  expect_equal(result[2, ], mad_score(m[2, ]))
})

test_that("robust_scale can scale across columns", {
  m <- matrix(c(1, 2, 2, 3, 100, 10, 20, 20, 30, 40), nrow = 2, byrow = TRUE)
  result <- robust_scale(m, margin = 2)

  expect_equal(result[, 1], mad_score(m[, 1]))
  expect_equal(result[, 5], mad_score(m[, 5]))
})

test_that("robust_scale can skip centering or scaling", {
  m <- matrix(c(1, 2, 2, 3, 100, 10, 20, 20, 30, 40), nrow = 2, byrow = TRUE)

  centred_only <- robust_scale(m, scale = FALSE)
  expect_equal(centred_only[1, ], m[1, ] - median(m[1, ]))

  scaled_only <- robust_scale(m, center = FALSE)
  expect_equal(scaled_only[1, ], m[1, ] / mad(m[1, ]))
})

test_that("robust_scale handles zero-MAD rows per `zero_mad`", {
  m <- matrix(c(5, 5, 5, 5, 1, 2, 2, 100), nrow = 2, byrow = TRUE)

  expect_equal(robust_scale(m, zero_mad = "zero")[1, ], rep(0, 4))
  expect_true(all(is.na(robust_scale(m, zero_mad = "na")[1, ])))
  expect_snapshot(robust_scale(m, zero_mad = "error"), error = TRUE)
})

test_that("robust_scale errors on non-matrix input", {
  expect_snapshot(robust_scale(data.frame(x = 1:5)), error = TRUE)
})

test_that("robust_scale errors on invalid margin", {
  m <- matrix(1:4, nrow = 2)
  expect_snapshot(robust_scale(m, margin = 3), error = TRUE)
})

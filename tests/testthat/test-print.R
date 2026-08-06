test_that("verbose summaries are concise and presentation-only", {
  d <- qc_data()
  expect_silent(mad_qc(d, c(lower = "lower", upper = "upper"), verbose = FALSE))
  expect_message(mad_qc(d, c(lower = "lower", upper = "upper"), verbose = TRUE), "observations checked")
  expect_message(mad_qc(d, c(lower = "lower", upper = "upper"), verbose = TRUE), "flagged by at least one metric")
  expect_message(mad_qc(d, c(lower = "lower", upper = "upper"), verbose = TRUE), "lower: 1 flagged")
})

test_that("verymad_qc printing returns invisibly and does not mutate", {
  r <- mad_qc(qc_data(), c(lower = "lower", upper = "upper"), output = "report", verbose = FALSE)
  before <- r
  printed <- capture.output(returned <- print(r))
  expect_match(paste(printed, collapse = "\n"), "<verymad_qc>")
  expect_match(paste(printed, collapse = "\n"), "Threshold status")
  expect_identical(returned, r)
  expect_identical(r, before)

  empty <- mad_qc(data.frame(a = numeric()), c(a = "upper"), output = "report", verbose = FALSE)
  expect_silent(capture.output(print(empty)))
})

test_that("the package exports exactly two functions", {
  expect_identical(sort(getNamespaceExports("veryMAD")), c("mad_qc", "mad_scale"))
})

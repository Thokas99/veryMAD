test_that("verbose summaries are concise and presentation-only", {
  d <- qc_data()
  expect_silent(mad_qc(d, c(lower = "lower", upper = "upper"), verbose = FALSE))
  expect_message(mad_qc(d, c(lower = "lower", upper = "upper"), verbose = TRUE), "observations checked")
  expect_message(mad_qc(d, c(lower = "lower", upper = "upper"), verbose = TRUE), "flagged by at least one metric")
  expect_message(mad_qc(d, c(lower = "lower", upper = "upper"), verbose = TRUE), "unresolved")
  expect_message(mad_qc(d, c(lower = "lower", upper = "upper"), verbose = TRUE), "without estimated thresholds")
  expect_message(mad_qc(data.frame(a = rep(NA_real_, 5)), c(a = "upper"),
                        verbose = TRUE), "5 unresolved")
})

test_that("verymad_qc printing returns invisibly and does not mutate", {
  r <- mad_qc(qc_data(), c(lower = "lower", upper = "upper"), output = "report", verbose = FALSE)
  before <- r
  printed <- capture.output(returned <- print(r))
  expect_match(paste(printed, collapse = "\n"), "<verymad_qc>")
  expect_match(paste(printed, collapse = "\n"), "Observations flagged")
  expect_match(paste(printed, collapse = "\n"), "Observations unresolved")
  expect_match(paste(printed, collapse = "\n"), "without estimated thresholds")
  expect_match(paste(printed, collapse = "\n"), "Threshold status")
  expect_identical(returned, r)
  expect_identical(r, before)

  empty <- mad_qc(data.frame(a = numeric()), c(a = "upper"), output = "report", verbose = FALSE)
  expect_silent(capture.output(print(empty)))

  unresolved <- mad_qc(
    data.frame(a = rep(NA_real_, 5)), c(a = "upper"),
    output = "report", verbose = FALSE
  )
  printed <- paste(capture.output(print(unresolved)), collapse = "\n")
  expect_match(printed, "Observations: 5")
  expect_match(printed, "Metrics: 1")
  expect_match(printed, "Observations flagged: 0")
  expect_match(printed, "Observations unresolved: 5")
  expect_match(printed, "Metrics without estimated thresholds: 1")
})

test_that("the package exports exactly two functions", {
  expect_identical(sort(getNamespaceExports("veryMAD")), c("mad_qc", "mad_scale"))
})

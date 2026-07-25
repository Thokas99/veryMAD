summary_fixture <- function() {
  data.frame(
    .obs = rep(1:3, 2), id = rep(c("a", "b", "c"), 2),
    sample = rep(c("S1", NA, "S1"), 2),
    metric = rep(c("x", "y"), each = 3), value = 1:6,
    lower = NA_real_, upper = 10, direction = rep(c("upper", "both"), each = 3),
    is_outlier = c(FALSE, TRUE, NA, TRUE, FALSE, FALSE),
    stringsAsFactors = FALSE
  )
}

test_that("metric summaries preserve first-seen groups and definitions", {
  qc <- summary_fixture()
  original <- qc
  out <- summarize_mad_qc(qc, "metric", group_by = "sample")
  expect_named(out, c("sample", "metric", "n_observations", "n_evaluated",
                      "n_missing", "n_outliers", "outlier_proportion", "direction"))
  expect_equal(out$sample, c("S1", NA, "S1", NA))
  expect_equal(out$metric, c("x", "x", "y", "y"))
  expect_equal(out$n_observations, c(2L, 1L, 2L, 1L))
  expect_equal(out$n_evaluated, c(1L, 1L, 2L, 1L))
  expect_equal(out$n_missing, c(1L, 0L, 0L, 0L))
  expect_equal(out$n_outliers, c(0L, 1L, 1L, 0L))
  expect_equal(out$outlier_proportion, c(0, 1, 0.5, 0))
  expect_s3_class(out, "tbl_df")
  expect_equal(qc, original)
  qc$batch <- rep(c("B1", "B2", "B1"), 2)
  grouped <- summarize_mad_qc(qc, "metric", group_by = c("sample", "batch"))
  expect_named(grouped, c("sample", "batch", "metric", "n_observations",
    "n_evaluated", "n_missing", "n_outliers", "outlier_proportion", "direction"))
})

test_that("observation summaries use three-valued overall flags", {
  out <- summarize_mad_qc(summary_fixture(), "observation", group_by = "sample")
  expect_named(out, c(".obs", "id", "sample", "n_metrics", "n_evaluated",
                      "n_missing", "n_outliers", "failed_metrics", "mad_qc_outlier"))
  expect_equal(out$.obs, 1:3)
  expect_equal(out$n_metrics, rep(2L, 3))
  expect_equal(out$failed_metrics, c("y", "x", ""))
  expect_equal(out$mad_qc_outlier, c(TRUE, TRUE, NA))
  expect_equal(out$n_missing, c(0L, 0L, 1L))
})

test_that("summary handles empty and malformed reports deliberately", {
  empty <- mad_qc(data.frame(x = numeric()), c(x = "both"))
  expect_equal(nrow(summarize_mad_qc(empty, "metric")), 0)
  expect_equal(nrow(summarize_mad_qc(empty, "observation")), 0)
  expect_s3_class(summarize_mad_qc(empty, "metric"), "tbl_df")
  expect_s3_class(summarize_mad_qc(empty, "observation"), "tbl_df")
  qc <- summary_fixture()
  expect_silent(summarize_mad_qc(as.data.frame(qc)))
  expect_error(summarize_mad_qc(qc[c(1, 1), ]), "duplicate")
  qc_bad_id <- qc
  qc_bad_id$id[4] <- "different"
  expect_error(summarize_mad_qc(qc_bad_id), "multiple identifiers")
  expect_error(summarize_mad_qc(qc[, setdiff(names(qc), "is_outlier")]), "must be")
  expect_error(summarize_mad_qc(qc, group_by = "batch"), "Missing grouping")
  all_missing <- qc[qc$metric == "x", ]
  all_missing$is_outlier <- NA
  expect_true(is.na(summarize_mad_qc(all_missing)$outlier_proportion))
  qc$direction[2] <- "lower"
  expect_error(summarize_mad_qc(qc), "inconsistent directions")
})

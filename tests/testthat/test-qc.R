test_that("mad_qc is tidy, grouped, transformed, and stable", {
  d <- data.frame(batch = c("a", NA, "a", NA), lane = c(1, 1, 1, 1),
                  counts = c(0, 9, 10, 99), mt = c(1, 2, 3, 40),
                  row.names = c("cell", "cell.1", "cell.2", "cell.3"))
  q <- mad_qc(d, c(counts = "both", mt = "upper"), nmads = 1,
              group_by = c("batch", "lane"), transform = c(counts = "log10p"))
  expect_equal(nrow(q), 8)
  expect_equal(q$.obs, rep(1:4, 2))
  expect_equal(q$batch[1:4], d$batch)
  expect_equal(q$value[1:4], log10(d$counts + 1))
  expect_equal(unique(q$direction[q$metric == "mt"]), "upper")
  expect_equal(length(unique(q$.obs[q$metric == "counts"])), 4)
})

test_that("mad_qc handles empty, singleton, and zero-MAD groups", {
  empty <- mad_qc(data.frame(x = numeric()), c(x = "both"))
  expect_equal(nrow(empty), 0)
  expect_s3_class(empty, "mad_qc")
  expect_named(annotate_mad_qc(data.frame(x = numeric()), empty),
               c("x", "x_mad_outlier", "mad_qc_outlier"))
  d <- data.frame(g = c("a", "b"), x = c(1, 2))
  expect_false(any(mad_qc(d, c(x = "both"), group_by = "g")$is_outlier))
  expect_true(all(is.na(mad_qc(d, c(x = "both"), group_by = "g", zero_mad = "na")$is_outlier)))
  expect_error(mad_qc(d, c(x = "both"), group_by = "g", zero_mad = "error"), "MAD is zero")
})

test_that("mad_qc validates metrics, groups, and transforms", {
  d <- data.frame(x = c(0, 1), text = c("a", "b"))
  expect_error(mad_qc(d, "both"), "named")
  expect_error(mad_qc(d, c(nope = "both")), "Missing metric")
  expect_error(mad_qc(d, c(text = "both")), "numeric")
  expect_error(mad_qc(d, c(x = "side")), "lower")
  expect_error(mad_qc(d, c(x = "both"), group_by = "nope"), "Missing grouping")
  expect_error(mad_qc(d, c(x = "both"), transform = c(nope = "log10")), "requested")
  expect_error(mad_qc(d, c(x = "both"), transform = c(x = "sqrt")), "Unsupported")
  expect_error(mad_qc(d, c(x = "both"), transform = c(x = "log10")), "positive")
  expect_equal(mad_qc(d, c(x = "both"), transform = c(x = "log10p"))$value, log10(d$x + 1))
})

test_that("annotation aligns by internal observation index", {
  d <- data.frame(x = c(1, 2, 100), y = c(1, 20, 3))
  q <- mad_qc(d, c(x = "upper", y = "upper"), nmads = 1)
  q <- q[rev(seq_len(nrow(q))), ]
  out <- annotate_mad_qc(d, q)
  expect_equal(out$x, d$x)
  expect_equal(out$x_mad_outlier, c(FALSE, FALSE, TRUE))
  expect_equal(out$mad_qc_outlier, out$x_mad_outlier | out$y_mad_outlier)
})

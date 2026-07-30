#' Bulk observation-level MAD quality control
#'
#' Rows are observations and `metrics` explicitly selects numeric QC columns.
#' A gene-by-sample expression matrix is not automatically a QC-metric table.
#' Selected metrics use `log1p` unless `transform` is `"none"` or supplies a
#' named partial override. Measurements are never changed or filtered.
#'
#' @param data A matrix, data frame, or tibble with observations in rows.
#' @param metrics Named character vector mapping columns to `"lower"`,
#'   `"upper"`, or `"both"`. Required; veryMAD never guesses metrics.
#' @param nmads Positive number of MADs used for thresholds.
#' @param transform `"log1p"`, `"none"`, `"identity"`, or `"log10"`, or a
#'   named partial override. Unnamed metrics in an override use `"log1p"`.
#' @param group_by Optional unique column names for within-group calculations.
#' @param verbose Print a concise QC summary?
#' @param overwrite Replace existing per-metric flag columns?
#' @return The input as a data frame with one logical `*_mad_outlier` column per
#'   metric. The compact report is available as `attr(result, "mad_qc")`.
#' @export
#' @examples
#' d <- data.frame(size = c(10, 11, 12, 100), rate = c(.9, .91, .92, .2))
#' out <- mad_qc_bulk(d, c(size = "upper", rate = "lower"),
#'                    transform = c(rate = "none"), verbose = FALSE)
#' attr(out, "mad_qc")$thresholds
mad_qc_bulk <- function(data, metrics = NULL, nmads = 3, transform = "log1p",
                        group_by = NULL, verbose = TRUE, overwrite = FALSE) {
  .arg_flag(verbose, "verbose"); .arg_flag(overwrite, "overwrite")
  if (!is.data.frame(data) && !is.matrix(data)) {
    cli::cli_abort(c("x" = "{.arg data} must be a matrix, data frame, or tibble."))
  }
  data <- as.data.frame(data, check.names = FALSE, stringsAsFactors = FALSE)
  metrics <- .validate_wrapper_metrics(data, metrics, "mad_qc_bulk")
  targets <- paste0(names(metrics), "_mad_outlier")
  conflicts <- intersect(targets, names(data))
  if (length(conflicts) && !overwrite) {
    cli::cli_abort(c("x" = "Flag column(s) already exist: {.field {conflicts}}.",
      "i" = "Use {.code overwrite = TRUE} to replace them."))
  }
  result <- .mad_qc_compact(data, metrics, nmads, transform, group_by)
  for (i in seq_along(targets)) data[[targets[[i]]]] <- result$flags[[names(metrics)[[i]]]]
  attr(data, "mad_qc") <- result
  if (verbose) .inform_mad_qc(result, "bulk QC", "observations")
  data
}

#' Single-cell observation-level MAD quality control
#'
#' Accepts cell metadata or a Seurat object. Annotation adds only
#' `mad_qc_outlier`; individual flags remain in the compact report. `TRUE` means
#' at least one metric failed, `FALSE` means all passed, and `NA` means none
#' failed but at least one could not be evaluated. Cells are never filtered.
#'
#' @inheritParams mad_qc_bulk
#' @param object A data frame, tibble, or Seurat object.
#' @param action Annotate the input or return the compact report.
#' @return For `action = "report"`, a compact `mad_qc_result`. Otherwise the
#'   metadata data frame or Seurat object with one `mad_qc_outlier` column. The
#'   data-frame annotation stores the report in `attr(result, "mad_qc")`.
#' @export
#' @examples
#' cells <- data.frame(counts = c(NA, 10, 11, 100), mt = c(2, 3, 4, 30))
#' mad_qc_sc(cells, c(counts = "both", mt = "upper"),
#'           transform = c(mt = "none"), action = "report", verbose = FALSE)
mad_qc_sc <- function(object, metrics = NULL, nmads = 3, transform = "log1p",
                      group_by = NULL, verbose = TRUE, overwrite = FALSE,
                      action = c("annotate", "report")) {
  action <- match.arg(action); .arg_flag(verbose, "verbose"); .arg_flag(overwrite, "overwrite")
  seurat <- inherits(object, "Seurat")
  if (seurat) {
    .require_namespace("SeuratObject", "`mad_qc_sc()` with a Seurat object")
    metadata <- object[[]]
  } else {
    if (!is.data.frame(object)) cli::cli_abort(c("x" = "{.arg object} must be cell metadata or a Seurat object."))
    metadata <- as.data.frame(object, check.names = FALSE, stringsAsFactors = FALSE)
  }
  if (action == "annotate" && "mad_qc_outlier" %in% names(metadata) && !overwrite) {
    cli::cli_abort(c("x" = "{.field mad_qc_outlier} already exists.",
      "i" = "Use {.code overwrite = TRUE} to replace it."))
  }
  metrics <- .validate_wrapper_metrics(metadata, metrics, "mad_qc_sc")
  result <- .mad_qc_compact(metadata, metrics, nmads, transform, group_by)
  if (verbose) .inform_mad_qc(result, "single-cell QC", "cells")
  if (action == "report") return(result)
  combined <- .overall_qc_flag(result$flags[names(metrics)])
  names(combined) <- rownames(metadata)
  if (seurat) return(SeuratObject::AddMetaData(object, metadata = combined, col.name = "mad_qc_outlier"))
  metadata$mad_qc_outlier <- combined
  attr(metadata, "mad_qc") <- result
  metadata
}

.validate_wrapper_metrics <- function(data, metrics, caller) {
  hint <- c("i" = "veryMAD will not guess which columns represent QC metrics.",
    ">" = "Provide a named vector mapping columns to {.val lower}, {.val upper}, or {.val both}.",
    "v" = "You choose the biology; veryMAD handles the MADness.")
  if (is.null(metrics)) cli::cli_abort(c("x" = "{.arg metrics} must be specified in {.fn {caller}}.", hint))
  if (!is.character(metrics) || !length(metrics) || is.null(names(metrics)) ||
      anyNA(names(metrics)) || any(names(metrics) == "") || anyDuplicated(names(metrics))) {
    cli::cli_abort(c("x" = "{.arg metrics} must be a non-empty named character vector with unique names.", hint))
  }
  if (!all(metrics %in% c("lower", "upper", "both"))) {
    cli::cli_abort(c("x" = "{.arg metrics} directions must be {.val lower}, {.val upper}, or {.val both}."))
  }
  missing <- setdiff(names(metrics), names(data))
  if (length(missing)) cli::cli_abort(c("x" = "Missing metric column{?s}: {.field {missing}}."))
  nonnumeric <- names(metrics)[!vapply(data[names(metrics)], is.numeric, logical(1))]
  if (length(nonnumeric)) cli::cli_abort(c("x" = "QC metric column{?s} must be numeric: {.field {nonnumeric}}."))
  metrics
}

.wrapper_transforms <- function(transform, metrics) {
  out <- stats::setNames(rep("log1p", length(metrics)), metrics)
  if (!is.character(transform) || !length(transform) || anyNA(transform)) {
    cli::cli_abort(c("x" = "{.arg transform} must be a supported scalar or named partial override."))
  }
  if (length(transform) == 1L && is.null(names(transform))) out[] <- transform else {
    if (is.null(names(transform)) || any(names(transform) == "") || anyDuplicated(names(transform))) {
      cli::cli_abort(c("x" = "A transformation vector must have unique metric names."))
    }
    unknown <- setdiff(names(transform), metrics)
    if (length(unknown)) cli::cli_abort(c("x" = "Transformation names must reference selected metrics: {.field {unknown}}."))
    out[names(transform)] <- transform
  }
  out[out == "identity"] <- "none"
  if (!all(out %in% c("log1p", "none", "log10"))) {
    cli::cli_abort(c("x" = "Transformations must be {.val log1p}, {.val none}, {.val identity}, or {.val log10}."))
  }
  out
}

.mad_qc_compact <- function(data, metrics, nmads, transform, group_by) {
  .arg_positive(nmads, "nmads")
  if (!is.null(group_by) && (!is.character(group_by) || !length(group_by) || anyNA(group_by) ||
      any(group_by == "") || anyDuplicated(group_by))) {
    cli::cli_abort(c("x" = "{.arg group_by} must be NULL or unique column names."))
  }
  missing_groups <- setdiff(group_by, names(data))
  if (length(missing_groups)) cli::cli_abort(c("x" = "Missing grouping column{?s}: {.field {missing_groups}}."))
  methods <- .wrapper_transforms(transform, names(metrics))
  n <- nrow(data); groups <- .group_index(data, group_by)
  flags <- data.frame(.obs = seq_len(n), id = .observation_ids(data), check.names = FALSE)
  threshold_rows <- list(); k <- 0L
  for (metric in names(metrics)) {
    raw <- data[[metric]]
    value <- .transform_metric(raw, if (methods[[metric]] == "none") "identity" else methods[[metric]], metric)
    metric_flags <- rep(NA, n)
    for (idx in split(seq_len(n), groups)) {
      stats <- .mad_stats(value[idx], 1.4826, TRUE, "zero")
      limits <- .mad_limits_from_stats(stats, nmads, metrics[[metric]], "zero")
      called <- .mad_flags(value[idx], stats, nmads, metrics[[metric]], "zero")
      metric_flags[idx] <- called; evaluated <- sum(!is.na(called)); k <- k + 1L
      row <- if (length(group_by)) data[idx[1L], group_by, drop = FALSE] else data.frame(.all = "all")
      row$metric <- metric; row$direction <- unname(metrics[[metric]]); row$transform <- methods[[metric]]
      row$n_observations <- length(idx); row$n_evaluated <- evaluated; row$n_missing <- sum(is.na(called))
      row$median <- stats$median; row$mad <- stats$mad; row$lower <- limits[["lower"]]; row$upper <- limits[["upper"]]
      inverse <- switch(methods[[metric]], log1p = expm1, log10 = function(x) 10^x, none = identity)
      row$lower_raw <- if (is.na(limits[["lower"]])) NA_real_ else inverse(limits[["lower"]])
      row$upper_raw <- if (is.na(limits[["upper"]])) NA_real_ else inverse(limits[["upper"]])
      row$n_outliers <- sum(called %in% TRUE)
      row$outlier_proportion <- if (evaluated) row$n_outliers / evaluated else NA_real_
      threshold_rows[[k]] <- row
    }
    flags[[metric]] <- metric_flags
  }
  thresholds <- if (length(threshold_rows)) do.call(rbind, threshold_rows) else {
    x <- data.frame(.all = character(), metric = character(), direction = character(), transform = character(),
      n_observations = integer(), n_evaluated = integer(), n_missing = integer(), median = numeric(), mad = numeric(),
      lower = numeric(), upper = numeric(), lower_raw = numeric(), upper_raw = numeric(), n_outliers = integer(),
      outlier_proportion = numeric(), stringsAsFactors = FALSE)
    if (length(group_by)) x <- cbind(data[FALSE, group_by, drop = FALSE], x[setdiff(names(x), ".all")])
    x
  }
  if (!length(group_by) && ".all" %in% names(thresholds)) thresholds$.all <- NULL
  out <- list(flags = tibble::as_tibble(flags), thresholds = tibble::as_tibble(thresholds),
    settings = list(metrics = metrics, transform = methods, nmads = nmads, group_by = group_by,
      observation_ids = .observation_ids(data)))
  class(out) <- c("mad_qc_result", "list")
  out
}

.qc_grade <- function(proportion) {
  if (is.na(proportion)) return("QC grade unavailable \u2014 no observations were evaluated.")
  if (proportion == 0) "Wonderful! :D" else if (proportion <= .05) "Looking good \u2014 only a few rebels." else
    if (proportion <= .15) "Not so good..." else "The wet lab is on fire! Metaphorically \u2014 inspect the QC."
}

.inform_mad_qc <- function(result, context, unit) {
  flags <- result$flags[names(result$settings$metrics)]
  overall <- .overall_qc_flag(flags); n <- nrow(flags); total <- sum(overall %in% TRUE)
  proportion <- if (n) total / n else NA_real_
  cli::cli_inform(c("{.strong veryMAD} \u2022 {context}",
    "v" = "{n} {unit} checked across {length(flags)} explicit metric{?s}"))
  shown <- utils::head(result$thresholds, 8L)
  for (i in seq_len(nrow(shown))) {
    x <- shown[i, ]; limit <- switch(x$direction,
      lower = paste0("< ", format(x$lower_raw, digits = 5)),
      upper = paste0("> ", format(x$upper_raw, digits = 5)),
      both = paste0("< ", format(x$lower_raw, digits = 5), " or > ", format(x$upper_raw, digits = 5)))
    cli::cli_inform(c("i" = "{x$metric} {x$direction} {limit} {x$transform} {x$n_outliers} flagged ({sprintf('%.1f%%', 100 * x$outlier_proportion)})"))
  }
  if (nrow(result$thresholds) > nrow(shown)) cli::cli_inform(c("i" = "Showing 8 threshold rows; see {.code result$thresholds} for all groups."))
  grade <- .qc_grade(proportion)
  cli::cli_inform(c("!" = "Overall: {total}/{n} {unit} flagged ({sprintf('%.1f%%', 100 * proportion)})",
    "v" = "{grade}"))
}

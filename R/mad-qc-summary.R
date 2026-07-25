#' Summarize an existing MAD QC report
#'
#' Summarizes flags already calculated by [mad_qc()]. Thresholds are not
#' recalculated.
#'
#' @param qc A tidy report created by [mad_qc()], either retaining its class or
#'   represented as an ordinary data frame with the required columns.
#' @param level Summarize by metric or original observation.
#' @param group_by `NULL` or one or more grouping column names in `qc`.
#' @return An ordinary data frame with deterministic columns and first-seen
#'   group, metric, or observation order.
#' @export
#' @examples
#' d <- data.frame(sample = c("a", "a", "b"), x = c(1, 2, 10))
#' qc <- mad_qc(d, c(x = "upper"), group_by = "sample")
#' summarize_mad_qc(qc, "metric", group_by = "sample")
#' summarize_mad_qc(qc, "observation")
summarize_mad_qc <- function(qc, level = c("metric", "observation"), group_by = NULL) {
  level <- match.arg(level)
  .validate_qc_report(qc, group_by)
  if (!nrow(qc)) return(.empty_qc_summary(level, qc, group_by))
  if (level == "metric") .summarize_qc_metrics(qc, group_by) else .summarize_qc_observations(qc, group_by)
}

.summarize_qc_metrics <- function(qc, group_by) {
  by <- c(group_by, "metric")
  group <- .group_index(qc, by)
  pieces <- lapply(split(seq_len(nrow(qc)), group), function(i) {
    flags <- qc$is_outlier[i]
    evaluated <- sum(!is.na(flags))
    out <- qc[i[1L], by, drop = FALSE]
    out$n_observations <- length(i)
    out$n_evaluated <- evaluated
    out$n_missing <- sum(is.na(flags))
    out$n_outliers <- sum(flags %in% TRUE)
    out$outlier_proportion <- if (evaluated) out$n_outliers / evaluated else NA_real_
    directions <- unique(qc$direction[i])
    if (length(directions) != 1L) stop("`qc` has inconsistent directions within a metric summary group.", call. = FALSE)
    out$direction <- directions
    out
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

.summarize_qc_observations <- function(qc, group_by) {
  group <- match(qc$.obs, unique(qc$.obs))
  pieces <- lapply(split(seq_len(nrow(qc)), group), function(i) {
    flags <- qc$is_outlier[i]
    failed <- qc$metric[i][flags %in% TRUE]
    out <- qc[i[1L], c(".obs", "id", group_by), drop = FALSE]
    if (length(group_by) && any(vapply(qc[i, group_by, drop = FALSE], function(x) length(unique(x)) > 1L, logical(1)))) {
      stop("`qc` has inconsistent grouping values within an observation.", call. = FALSE)
    }
    out$n_metrics <- length(i)
    out$n_evaluated <- sum(!is.na(flags))
    out$n_missing <- sum(is.na(flags))
    out$n_outliers <- length(failed)
    out$failed_metrics <- paste(failed, collapse = ", ")
    out$mad_qc_outlier <- length(failed) > 0L
    out
  })
  out <- do.call(rbind, pieces)
  rownames(out) <- NULL
  out
}

.empty_qc_summary <- function(level, qc, group_by) {
  groups <- qc[group_by]
  if (level == "metric") {
    return(cbind(groups, data.frame(metric = character(), n_observations = integer(),
      n_evaluated = integer(), n_missing = integer(), n_outliers = integer(),
      outlier_proportion = numeric(), direction = character(), stringsAsFactors = FALSE)))
  }
  cbind(data.frame(.obs = numeric(), id = character(), stringsAsFactors = FALSE), groups,
    data.frame(n_metrics = integer(), n_evaluated = integer(), n_missing = integer(),
      n_outliers = integer(), failed_metrics = character(), mad_qc_outlier = logical(),
      stringsAsFactors = FALSE))
}

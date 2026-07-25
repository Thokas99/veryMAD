#' Add MAD QC flags to a data frame
#'
#' @param data Original data frame used by [mad_qc()].
#' @param qc A tidy result from [mad_qc()].
#' @param suffix Suffix for per-metric logical columns.
#' @param overall Name of the overall logical flag.
#' @param overwrite Replace only conflicting veryMAD flag columns?
#' @return `data` with aligned logical QC columns appended.
#' @export
#' @examples
#' d <- data.frame(x = c(1, 1, 1, 10))
#' q <- mad_qc(d, c(x = "upper"))
#' annotate_mad_qc(d, q)
annotate_mad_qc <- function(data, qc, suffix = "_mad_outlier", overall = "mad_qc_outlier",
                            overwrite = FALSE) {
  if (!is.data.frame(data) || !is.data.frame(qc) || !all(c(".obs", "metric", "is_outlier") %in% names(qc))) {
    stop("`data` must be a data frame and `qc` a result from `mad_qc()`.", call. = FALSE)
  }
  if (!is.character(suffix) || length(suffix) != 1L || is.na(suffix) ||
      !is.character(overall) || length(overall) != 1L || is.na(overall) || overall == "") {
    stop("`suffix` and `overall` must be single nonmissing strings.", call. = FALSE)
  }
  .arg_flag(overwrite, "overwrite")
  .validate_qc_report(qc)
  metrics <- if (nrow(qc)) unique(qc$metric) else attr(qc, "metrics")
  flags <- stats::setNames(vector("list", length(metrics)), paste0(metrics, suffix))
  intended <- c(names(flags), overall)
  if (anyDuplicated(intended)) stop("Per-metric and overall flag column names must be unique.", call. = FALSE)
  conflicts <- intersect(intended, names(data))
  if (length(conflicts) && !overwrite) {
    stop(sprintf("Flag column(s) already exist: %s. Use `overwrite = TRUE` to replace them.",
      paste(conflicts, collapse = ", ")), call. = FALSE)
  }
  for (i in seq_along(metrics)) {
    part <- qc[qc$metric == metrics[[i]], c(".obs", "id", "is_outlier")]
    if (nrow(part) != nrow(data) || !setequal(part$.obs, seq_len(nrow(data)))) stop("`qc` does not align with `data`.", call. = FALSE)
    order <- match(seq_len(nrow(data)), part$.obs)
    if (!identical(as.character(part$id[order]), .observation_ids(data))) stop("`qc` identifiers do not align with `data`.", call. = FALSE)
    flags[[i]] <- part$is_outlier[order]
  }
  for (name in names(flags)) data[[name]] <- flags[[name]]
  data[[overall]] <- if (length(flags)) .overall_qc_flag(as.data.frame(flags)) else rep(FALSE, nrow(data))
  data
}

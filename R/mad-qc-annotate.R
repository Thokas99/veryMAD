#' Add MAD QC flags to a data frame
#'
#' @param data Original data frame used by [mad_qc()].
#' @param qc A tidy result from [mad_qc()].
#' @param suffix Suffix for per-metric logical columns.
#' @param overall Name of the overall logical flag.
#' @return `data` with aligned logical QC columns appended.
#' @export
#' @examples
#' d <- data.frame(x = c(1, 1, 1, 10))
#' q <- mad_qc(d, c(x = "upper"))
#' annotate_mad_qc(d, q)
annotate_mad_qc <- function(data, qc, suffix = "_mad_outlier", overall = "mad_qc_outlier") {
  if (!is.data.frame(data) || !is.data.frame(qc) || !all(c(".obs", "metric", "is_outlier") %in% names(qc))) {
    stop("`data` must be a data frame and `qc` a result from `mad_qc()`.", call. = FALSE)
  }
  if (!is.character(suffix) || length(suffix) != 1L || is.na(suffix) ||
      !is.character(overall) || length(overall) != 1L || is.na(overall) || overall == "") {
    stop("`suffix` and `overall` must be single nonmissing strings.", call. = FALSE)
  }
  metrics <- if (nrow(qc)) unique(qc$metric) else attr(qc, "metrics")
  flags <- stats::setNames(vector("list", length(metrics)), paste0(metrics, suffix))
  for (i in seq_along(metrics)) {
    part <- qc[qc$metric == metrics[[i]], c(".obs", "is_outlier")]
    if (nrow(part) != nrow(data) || !setequal(part$.obs, seq_len(nrow(data)))) stop("`qc` does not align with `data`.", call. = FALSE)
    flags[[i]] <- part$is_outlier[match(seq_len(nrow(data)), part$.obs)]
  }
  for (name in names(flags)) data[[name]] <- flags[[name]]
  data[[overall]] <- if (length(flags)) Reduce(function(x, y) x | y, flags) else rep(FALSE, nrow(data))
  data
}

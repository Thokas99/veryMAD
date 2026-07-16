#' Flag QC outliers using directional MAD thresholds
#'
#' Applies a MAD-based threshold to one or more QC metrics (e.g. read
#' counts, feature counts, percent mitochondrial reads), optionally within
#' groups (samples, batches), and returns a tidy audit table rather than
#' silently dropping observations.
#'
#' @param data A data frame, one row per observation (e.g. cell or sample).
#' @param metrics A named character vector. Names are columns in `data`;
#'   values give the outlier direction for that metric, one of `"lower"`,
#'   `"upper"`, or `"both"`.
#' @param nmads A positive number of MADs from the median beyond which an
#'   observation is flagged.
#' @param transform An optional named character vector giving a transform
#'   to apply to specific metrics before computing thresholds, one of
#'   `"identity"` (the default for any metric not listed) or `"log10"`.
#' @param group_by An optional column name in `data` (e.g. sample or batch)
#'   within which medians and MADs are computed separately.
#'
#' @return A data frame with one row per observation-metric combination,
#'   with columns `id` (the row name or index of `data`), `metric`,
#'   `value` (possibly transformed), `median`, `mad`, `lower`, `upper`, and
#'   `is_outlier`.
#' @export
#'
#' @examples
#' cell_metadata <- data.frame(
#'   nCount_RNA = c(500, 600, 550, 20000, 580),
#'   percent.mt = c(2, 3, 2.5, 3, 40),
#'   sample = c("a", "a", "a", "a", "a")
#' )
#' mad_qc(
#'   cell_metadata,
#'   metrics = c(nCount_RNA = "lower", percent.mt = "upper")
#' )
mad_qc <- function(
  data,
  metrics,
  nmads = 3,
  transform = NULL,
  group_by = NULL
) {
  if (!is.data.frame(data)) {
    cli::cli_abort("{.arg data} must be a data frame.")
  }

  metric_names <- names(metrics)
  if (is.null(metric_names) || any(metric_names == "")) {
    cli::cli_abort("{.arg metrics} must be a named character vector.")
  }
  if (!all(metrics %in% c("lower", "upper", "both"))) {
    cli::cli_abort(
      '{.arg metrics} values must be "lower", "upper", or "both".'
    )
  }
  missing_cols <- setdiff(metric_names, names(data))
  if (length(missing_cols) > 0) {
    cli::cli_abort("Column(s) {.val {missing_cols}} not found in {.arg data}.")
  }

  if (!is.numeric(nmads) || length(nmads) != 1L || nmads <= 0) {
    cli::cli_abort("{.arg nmads} must be one positive number.")
  }

  if (!is.null(transform)) {
    if (!all(transform %in% c("identity", "log10"))) {
      cli::cli_abort(
        '{.arg transform} values must be "identity" or "log10".'
      )
    }
  }

  if (!is.null(group_by)) {
    if (!group_by %in% names(data)) {
      cli::cli_abort("Column {.val {group_by}} not found in {.arg data}.")
    }
    groups <- data[[group_by]]
  } else {
    groups <- rep(1L, nrow(data))
  }

  id <- if (!is.null(rownames(data))) rownames(data) else seq_len(nrow(data))

  results <- lapply(metric_names, function(metric) {
    direction <- metrics[[metric]]
    transform_type <- if (metric %in% names(transform)) {
      transform[[metric]]
    } else {
      "identity"
    }

    raw <- data[[metric]]
    value <- switch(
      transform_type,
      log10 = log10(raw),
      identity = raw
    )

    thresholds <- lapply(split(seq_along(value), groups), function(idx) {
      med <- stats::median(value[idx], na.rm = TRUE)
      m <- stats::mad(value[idx], center = med, na.rm = TRUE)
      lower <- if (direction %in% c("lower", "both")) med - nmads * m else -Inf
      upper <- if (direction %in% c("upper", "both")) med + nmads * m else Inf
      data.frame(idx = idx, median = med, mad = m, lower = lower, upper = upper)
    })
    thresholds <- do.call(rbind, thresholds)
    thresholds <- thresholds[order(thresholds$idx), ]

    data.frame(
      id = id,
      metric = metric,
      value = value,
      median = thresholds$median,
      mad = thresholds$mad,
      lower = thresholds$lower,
      upper = thresholds$upper,
      is_outlier = value < thresholds$lower | value > thresholds$upper
    )
  })

  out <- do.call(rbind, results)
  rownames(out) <- NULL
  out
}

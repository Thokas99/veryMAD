#' Flag observations with explicit MAD thresholds
#'
#' `mad_qc()` calculates median absolute deviation (MAD) thresholds for selected
#' numeric metadata columns and returns a long, auditable report. Rows in `data`
#' are observations: bulk RNA-seq libraries, single cells, spatial spots, or any
#' other assay-level unit with observation metadata.
#'
#' Count-like QC metrics such as library size, total RNA counts, and detected
#' features are often strongly right-skewed. In those cases, estimating MAD
#' thresholds after an explicit log transformation such as `log10p` can give a
#' more useful reference distribution. Bounded percentages, proportions, and rates
#' are usually kept on their original scale because a log transform changes the
#' interpretation of already bounded measurements.
#'
#' Transformations are explicit and per metric. `transform = NULL` is equivalent
#' to identity transformation for every metric, and `mad_qc()` never infers a
#' transformation from a metric name. Thresholds, medians, MADs, and the `value`
#' column are expressed on the calculation scale. The `raw_value` column always
#' preserves the exact input measurement scale.
#'
#' `direction = "both"` flags both tails of the selected metric, but it does
#' not assign a biological cause. For example, an upper-tail count or feature flag
#' means unusually high relative to the reference distribution; it is not a
#' doublet call or an automatic filtering decision.
#'
#' @param data A data frame containing one observation per row.
#' @param metrics A named character vector mapping metric names in `data` to one
#'   of `"lower"`, `"upper"`, or `"both"`.
#' @param nmads Number of MADs from the median used to define thresholds.
#' @param group_by Optional column name used to estimate thresholds within groups.
#' @param transform Optional named character vector mapping metrics to explicit
#'   transformations. Supported values are `"identity"`, `"log10"`, and
#'   `"log10p"`. Metrics absent from this vector use identity transformation.
#' @param constant Consistency constant passed to MAD calculation.
#' @param na_rm Logical; remove missing values before threshold estimation.
#' @param zero_mad How to handle groups with zero MAD: `"zero"`, `"na"`, or
#'   `"error"`.
#'
#' @return A `mad_qc` data frame with one row per observation-metric pair. It
#'   includes `raw_value`, transformed `value`, threshold columns, `direction`, and
#'   `is_outlier`.
#' @export
#'
#' @examples
#' bulk_metadata <- data.frame(
#'   sample_id = paste0("sample_", 1:8),
#'   library_size = c(2e6, 3e7, 3.2e7, 3.4e7, 3.1e7, 3.3e7, 3.5e7, 3.6e7),
#'   mapping_rate = c(0.55, 0.92, 0.94, 0.93, 0.95, 0.94, 0.93, 0.92)
#' )
#' mad_qc(
#'   bulk_metadata,
#'   metrics = c(library_size = "lower", mapping_rate = "lower"),
#'   transform = c(library_size = "log10p")
#' )
#' cell_metadata <- data.frame(
#'   cell_id = paste0("cell_", 1:8),
#'   nCount_RNA = c(200, 8000, 8500, 9000, 8700, 9200, 9500, 70000),
#'   nFeature_RNA = c(150, 2800, 3000, 3100, 2950, 3200, 3300, 8000),
#'   percent.mt = c(3, 4, 5, 4, 6, 5, 4, 18)
#' )
#' mad_qc(
#'   cell_metadata,
#'   metrics = c(nCount_RNA = "both", nFeature_RNA = "both", percent.mt = "upper"),
#'   transform = c(nCount_RNA = "log10p", nFeature_RNA = "log10p")
#' )
mad_qc <- function(data, metrics, nmads = 3, group_by = NULL, transform = NULL,
                   constant = 1.4826, na_rm = TRUE,
                   zero_mad = c("zero", "na", "error")) {
  if (!is.data.frame(data)) stop("`data` must be a data frame.", call. = FALSE)
  zero_mad <- match.arg(zero_mad)
  .arg_positive(nmads, "nmads"); .arg_positive(constant, "constant"); .arg_flag(na_rm, "na_rm")
  if (!is.character(metrics) || !length(metrics) || is.null(names(metrics)) ||
      anyNA(names(metrics)) || any(names(metrics) == "") || anyDuplicated(names(metrics))) {
    stop("`metrics` must be a nonempty named character vector with unique names.", call. = FALSE)
  }
  if (!all(metrics %in% c("lower", "upper", "both"))) {
    stop('`metrics` values must be "lower", "upper", or "both".', call. = FALSE)
  }
  metric_names <- names(metrics)
  missing_metrics <- setdiff(metric_names, names(data))
  if (length(missing_metrics)) stop(sprintf("Missing metric column(s): %s.", paste(missing_metrics, collapse = ", ")), call. = FALSE)
  if (!all(vapply(data[metric_names], is.numeric, logical(1)))) stop("All QC metric columns must be numeric.", call. = FALSE)
  if (!is.null(group_by) && (!is.character(group_by) || !length(group_by) || anyNA(group_by) ||
      any(group_by == "") || anyDuplicated(group_by))) {
    stop("`group_by` must be NULL or one or more unique column names.", call. = FALSE)
  }
  missing_groups <- setdiff(group_by, names(data))
  if (length(missing_groups)) stop(sprintf("Missing grouping column(s): %s.", paste(missing_groups, collapse = ", ")), call. = FALSE)
  transform <- .validate_transforms(transform, metric_names)
  n <- nrow(data)
  if (!n) return(.empty_qc(data, group_by, metric_names, transform, nmads, constant, zero_mad))
  group <- .group_index(data, group_by)
  ids <- .observation_ids(data)
  pieces <- lapply(metric_names, function(metric) {
    raw <- data[[metric]]
    value <- .transform_metric(raw, transform[[metric]], metric)
    median <- mad <- lower <- upper <- rep(NA_real_, n)
    is_outlier <- rep(NA, n)
    for (i in split(seq_len(n), group)) {
      stats <- .mad_stats(value[i], constant, na_rm, zero_mad)
      limits <- .mad_limits_from_stats(stats, nmads, metrics[[metric]], zero_mad)
      median[i] <- stats$median
      mad[i] <- stats$mad
      lower[i] <- limits[["lower"]]
      upper[i] <- limits[["upper"]]
      is_outlier[i] <- .mad_flags(value[i], stats, nmads, metrics[[metric]], zero_mad)
    }
    out <- data.frame(.obs = seq_len(n), id = ids, metric = metric,
                      raw_value = raw, value = value,
                      median = median, mad = mad, lower = lower, upper = upper,
                      direction = unname(metrics[[metric]]),
                      is_outlier = is_outlier,
                      stringsAsFactors = FALSE, check.names = FALSE)
    if (length(group_by)) out <- cbind(out[c(".obs", "id")], data[group_by], out[setdiff(names(out), c(".obs", "id"))])
    out
  })
  out <- .as_mad_qc(do.call(rbind, pieces))
  attr(out, "metrics") <- metric_names
  attr(out, "group_by") <- group_by
  attr(out, "transform") <- transform
  attr(out, "nmads") <- nmads
  attr(out, "constant") <- constant
  attr(out, "zero_mad") <- zero_mad
  attr(out, "observation_metadata") <- data[setdiff(names(data), metric_names)]
  out
}

.observation_ids <- function(data) {
  if (.row_names_info(data, type = 1L) < 0L) as.character(seq_len(nrow(data))) else rownames(data)
}

.validate_transforms <- function(transform, metrics) {
  out <- stats::setNames(rep("identity", length(metrics)), metrics)
  if (is.null(transform)) return(out)
  if (!is.character(transform) || !length(transform) || is.null(names(transform)) ||
      any(names(transform) == "") || anyDuplicated(names(transform))) {
    stop("`transform` must be a named character vector.", call. = FALSE)
  }
  if (!all(names(transform) %in% metrics)) stop("`transform` names must reference requested metrics.", call. = FALSE)
  if (!all(transform %in% c("identity", "log10", "log10p"))) stop("Unsupported transformation.", call. = FALSE)
  out[names(transform)] <- transform
  out
}

.transform_metric <- function(x, method, metric) {
  finite <- x[!is.na(x)]
  if (any(!is.finite(finite))) stop(sprintf("Metric `%s` contains non-finite values.", metric), call. = FALSE)
  if (method == "log10" && any(finite <= 0)) stop(sprintf("`log10` requires positive values in `%s`.", metric), call. = FALSE)
  if (method == "log10p" && any(finite < 0)) stop(sprintf("`log10p` requires nonnegative values in `%s`.", metric), call. = FALSE)
  switch(method, identity = x, log10 = log10(x), log10p = log10(x + 1))
}

.group_index <- function(data, group_by) {
  if (!length(group_by)) return(rep.int(1L, nrow(data)))
  keys <- lapply(data[group_by], function(x) {
    x <- as.character(x)
    ifelse(is.na(x), "M", paste0("V", nchar(x), ":", x))
  })
  match(do.call(paste, c(keys, sep = "\r")), unique(do.call(paste, c(keys, sep = "\r"))))
}

.empty_qc <- function(data, group_by, metrics, transform, nmads, constant, zero_mad) {
  out <- data.frame(.obs = integer(), id = character(), stringsAsFactors = FALSE)
  if (length(group_by)) out <- cbind(out, data[group_by])
  out <- cbind(out, data.frame(metric = character(), raw_value = numeric(), value = numeric(),
    median = numeric(), mad = numeric(), lower = numeric(), upper = numeric(),
    direction = character(), is_outlier = logical(), stringsAsFactors = FALSE))
  out <- .as_mad_qc(out)
  attr(out, "metrics") <- metrics
  attr(out, "group_by") <- group_by
  attr(out, "transform") <- transform
  attr(out, "nmads") <- nmads
  attr(out, "constant") <- constant
  attr(out, "zero_mad") <- zero_mad
  attr(out, "observation_metadata") <- data[setdiff(names(data), metrics)]
  out
}

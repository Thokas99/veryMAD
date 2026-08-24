#' Explicit MAD quality-control annotation
#'
#' Rows are observations and `metrics` explicitly selects numeric QC columns.
#' One reference distribution is calculated per metric across all observations;
#' observations are never filtered automatically.
#'
#' @param data A data frame, tibble, numeric matrix, or Seurat object.
#' @param metrics A non-empty named character vector mapping columns to
#'   `"lower"`, `"upper"`, or `"both"`.
#' @param nmads Positive number of MADs from the median for the limits.
#' @param transform A scalar transformation or named partial overrides.
#'   Supported values are `"none"`, `"log1p"`, and `"log10"`.
#' @param output Return annotated data (`"annotate"`) or a compact report
#'   (`"report"`).
#' @param min_n Minimum number of finite, non-missing observations required.
#' @param zero_mad Handling of a zero MAD: `"na"`, `"zero"`, or `"error"`.
#' @param overwrite Replace existing generated flag columns?
#' @param verbose Print a concise neutral summary?
#' @return Annotated input or a `verymad_qc` list with `flags`, `thresholds`,
#'   and `settings`.
#' @export
#' @examples
#' qc_metadata <- data.frame(
#'   sample = paste0("sample_", 1:6),
#'   library_size = c(2.4e6, 2.5e6, 2.6e6, 2.7e6, 0.5e6, 2.5e6),
#'   pct_mito = c(.04, .05, .06, .05, .07, .30)
#' )
#' qc <- mad_qc(
#'   qc_metadata,
#'   metrics = c(library_size = "lower", pct_mito = "upper"),
#'   verbose = FALSE
#' )
#' qc[, c("sample", "mad_qc_outlier")]
mad_qc <- function(data, metrics, nmads = 3, transform = "none",
                   output = c("annotate", "report"), min_n = 5,
                   zero_mad = c("na", "zero", "error"),
                   overwrite = FALSE, verbose = TRUE) {
  output <- match.arg(output)
  zero_mad <- match.arg(zero_mad)
  .validate_scalar_flag(overwrite, "overwrite")
  .validate_scalar_flag(verbose, "verbose")
  .validate_positive(nmads, "nmads")
  .validate_positive(min_n, "min_n")
  if (min_n != as.integer(min_n)) stop("`min_n` must be a positive integer.", call. = FALSE)

  is_seurat <- inherits(data, "Seurat")
  original <- data
  if (is_seurat) {
    if (!requireNamespace("SeuratObject", quietly = TRUE)) {
      stop("Seurat input requires the `SeuratObject` package.", call. = FALSE)
    }
    data <- data[[]]
  } else if (is.matrix(data)) {
    if (!is.numeric(data)) stop("`data` must be a numeric matrix.", call. = FALSE)
    data <- as.data.frame(data, check.names = FALSE, stringsAsFactors = FALSE)
  } else if (!is.data.frame(data)) {
    stop("`data` must be a data frame, numeric matrix, or Seurat object.", call. = FALSE)
  }

  metrics <- .validate_metrics(data, metrics)
  transforms <- .resolve_transforms(transform, names(metrics))
  targets <- c(paste0(names(metrics), "_mad_outlier"), "mad_qc_outlier")
  if (output == "annotate" && !overwrite) {
    conflicts <- intersect(targets, names(data))
    if (length(conflicts)) stop(sprintf("Flag column(s) already exist: %s. Use `overwrite = TRUE`.",
      paste(conflicts, collapse = ", ")), call. = FALSE)
  }

  result <- .qc_engine(data, metrics, transforms, nmads, min_n, zero_mad)
  annotated <- result$annotated
  result$annotated <- NULL
  if (verbose) .inform_qc(result, nrow(data))
  if (output == "report") return(result)

  if (is_seurat) {
    metadata <- result$flags[setdiff(names(result$flags), "id")]
    names(metadata)[seq_along(metrics)] <- paste0(names(metrics), "_mad_outlier")
    rownames(metadata) <- result$flags$id
    return(SeuratObject::AddMetaData(original, metadata = metadata))
  }
  if (is.matrix(original)) return(annotated)
  out <- data
  for (nm in names(annotated)) if (nm %in% targets) out[[nm]] <- annotated[[nm]]
  out
}

.validate_scalar_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
}

.validate_positive <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0) stop(sprintf("`%s` must be positive and finite.", name), call. = FALSE)
}

.validate_metrics <- function(data, metrics) {
  if (!is.character(metrics) || !length(metrics) || is.null(names(metrics)) ||
      anyNA(names(metrics)) || any(names(metrics) == "") || anyDuplicated(names(metrics))) {
    stop("`metrics` must be a non-empty named character vector with unique names.", call. = FALSE)
  }
  if (!all(metrics %in% c("lower", "upper", "both"))) stop("Metric directions must be `lower`, `upper`, or `both`.", call. = FALSE)
  missing <- setdiff(names(metrics), names(data))
  if (length(missing)) stop(sprintf("Missing metric column(s): %s.", paste(missing, collapse = ", ")), call. = FALSE)
  bad <- names(metrics)[!vapply(data[names(metrics)], is.numeric, logical(1))]
  if (length(bad)) stop(sprintf("QC metric column(s) must be numeric: %s.", paste(bad, collapse = ", ")), call. = FALSE)
  invisible(metrics)
}

.resolve_transforms <- function(transform, metrics) {
  out <- stats::setNames(rep("none", length(metrics)), metrics)
  if (!is.character(transform) || !length(transform) || anyNA(transform)) stop("`transform` must be a supported scalar or named partial override.", call. = FALSE)
  if (length(transform) == 1L && is.null(names(transform))) {
    out[] <- transform
  } else {
    if (is.null(names(transform)) || any(names(transform) == "") || anyDuplicated(names(transform))) stop("Named transformations must have unique metric names.", call. = FALSE)
    unknown <- setdiff(names(transform), metrics)
    if (length(unknown)) stop(sprintf("Transformation names must reference selected metrics: %s.", paste(unknown, collapse = ", ")), call. = FALSE)
    out[names(transform)] <- transform
  }
  out[out == "identity"] <- "none"
  if (!all(out %in% c("none", "log1p", "log10"))) stop("Transformations must be `none`, `log1p`, or `log10`.", call. = FALSE)
  out
}

.transform_metric <- function(x, method, metric) {
  if (any(is.nan(x) | is.infinite(x))) stop(sprintf("Metric `%s` contains non-finite values.", metric), call. = FALSE)
  finite <- x[!is.na(x)]
  if (method == "log1p" && any(finite < 0)) stop(sprintf("`log1p` requires non-negative values in `%s`.", metric), call. = FALSE)
  if (method == "log10" && any(finite <= 0)) stop(sprintf("`log10` requires positive values in `%s`.", metric), call. = FALSE)
  switch(method, none = x, log1p = log1p(x), log10 = log10(x))
}

.observation_ids <- function(data) {
  ids <- rownames(data)
  if (is.null(ids) || length(ids) != nrow(data) || anyNA(ids) || any(ids == "") || anyDuplicated(ids)) as.character(seq_len(nrow(data))) else ids
}

.qc_engine <- function(data, metrics, transforms, nmads, min_n, zero_mad) {
  n <- nrow(data); ids <- .observation_ids(data)
  flags <- data.frame(id = ids, stringsAsFactors = FALSE, check.names = FALSE)
  thresholds <- vector("list", length(metrics)); insufficient <- character()
  for (i in seq_along(metrics)) {
    metric <- names(metrics)[i]; direction <- unname(metrics[[i]])
    raw <- data[[metric]]; value <- .transform_metric(raw, transforms[[metric]], metric)
    usable <- value[is.finite(value) & !is.na(value)]
    status <- "ok"; centre <- spread <- lower <- upper <- NA_real_
    if (!length(usable)) status <- "all_missing" else if (length(usable) < min_n) {
      status <- "insufficient_n"; insufficient <- c(insufficient, metric)
    } else {
      centre <- stats::median(usable); spread <- stats::mad(usable, center = centre, constant = 1.4826)
      if (spread == 0) {
        status <- "zero_mad"
        if (zero_mad == "error") stop(sprintf("MAD is zero for metric `%s`.", metric), call. = FALSE)
        if (zero_mad == "zero") {
          if (direction %in% c("lower", "both")) lower <- centre
          if (direction %in% c("upper", "both")) upper <- centre
        }
      } else {
        if (direction %in% c("lower", "both")) lower <- centre - nmads * spread
        if (direction %in% c("upper", "both")) upper <- centre + nmads * spread
      }
    }
    flag <- rep(NA, n)
    if (status == "ok" || (status == "zero_mad" && zero_mad == "zero")) {
      present <- !is.na(value)
      flag[present] <- switch(direction, lower = value[present] < lower, upper = value[present] > upper, both = value[present] < lower | value[present] > upper)
    }
    flags[[metric]] <- flag
    inverse <- switch(transforms[[metric]], none = identity, log1p = expm1, log10 = function(x) 10^x)
    thresholds[[i]] <- data.frame(metric = metric, direction = direction, transform = transforms[[metric]],
      median = centre, mad = spread, lower = lower, upper = upper,
      lower_raw = if (is.na(lower)) NA_real_ else inverse(lower),
      upper_raw = if (is.na(upper)) NA_real_ else inverse(upper), status = status,
      stringsAsFactors = FALSE, check.names = FALSE)
  }
  if (length(insufficient)) warning(sprintf("Insufficient usable observations for metric(s): %s.", paste(insufficient, collapse = ", ")), call. = FALSE)
  flags$mad_qc_outlier <- .combine_flags(flags[names(metrics)])
  annotated <- data
  for (metric in names(metrics)) annotated[[paste0(metric, "_mad_outlier")]] <- flags[[metric]]
  annotated$mad_qc_outlier <- flags$mad_qc_outlier
  out <- list(flags = flags, thresholds = do.call(rbind, thresholds), settings = list(
    metrics = metrics, transform = transforms, nmads = nmads, constant = 1.4826,
    min_n = min_n, zero_mad = zero_mad))
  class(out) <- c("verymad_qc", "list")
  out$annotated <- annotated
  out
}

.combine_flags <- function(flags) {
  if (!nrow(flags)) return(logical())
  vapply(seq_len(nrow(flags)), function(i) {
    x <- unlist(flags[i, ], use.names = FALSE)
    if (any(x %in% TRUE)) TRUE else if (anyNA(x)) NA else FALSE
  }, logical(1))
}

.inform_qc <- function(result, n) {
  metrics <- names(result$settings$metrics); combined <- result$flags$mad_qc_outlier
  cli::cli_inform("veryMAD \u2022 observation QC")
  cli::cli_inform(sprintf("%d observations checked across %d metrics", n, length(metrics)))
  cli::cli_inform(sprintf("%d observations flagged by at least one metric", sum(combined %in% TRUE)))
  for (metric in metrics) {
    cli::cli_inform(sprintf("%s: %d flagged", metric, sum(result$flags[[metric]] %in% TRUE)))
  }
}

#' @method print verymad_qc
#' @export
print.verymad_qc <- function(x, ...) {
  metrics <- setdiff(names(x$flags), c("id", "mad_qc_outlier"))
  cat("<verymad_qc>\n")
  cat("Observations:", nrow(x$flags), "\n")
  cat("Metrics:", length(metrics), "\n")
  cat("Flagged by any metric:", sum(x$flags$mad_qc_outlier %in% TRUE), "\n")
  if (length(metrics)) cat("Metric names:", paste(metrics, collapse = ", "), "\n")
  cat("\nThreshold status:\n")
  statuses <- table(x$thresholds$status)
  for (status in names(statuses)) cat("  ", status, ": ", statuses[[status]], "\n", sep = "")
  invisible(x)
}

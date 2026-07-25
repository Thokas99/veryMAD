.arg_positive <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0) {
    stop(sprintf("`%s` must be one finite positive number.", name), call. = FALSE)
  }
}

.arg_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
}

.numeric_vector <- function(x) {
  if (!is.numeric(x) || !is.null(dim(x))) {
    stop("`x` must be a numeric vector.", call. = FALSE)
  }
  if (any(is.infinite(x), na.rm = TRUE)) {
    stop("`x` must not contain Inf or -Inf.", call. = FALSE)
  }
  if (any(is.nan(x))) stop("`x` must not contain NaN.", call. = FALSE)
}

.numeric_matrix <- function(x, allow_sparse = FALSE) {
  if (inherits(x, "sparseMatrix")) {
    if (!allow_sparse) {
      stop("Sparse matrices are not supported here because median centering can destroy sparsity; convert explicitly only if safe.", call. = FALSE)
    }
    if (!inherits(x, "dgCMatrix")) {
      stop("Only Matrix::dgCMatrix sparse matrices are supported.", call. = FALSE)
    }
    values <- x@x
    if (any(is.infinite(values), na.rm = TRUE)) stop("`x` must not contain Inf or -Inf.", call. = FALSE)
    if (any(is.nan(values))) stop("`x` must not contain NaN.", call. = FALSE)
    return(invisible(TRUE))
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`x` must be an ordinary numeric matrix.", call. = FALSE)
  }
  if (any(is.infinite(x), na.rm = TRUE)) {
    stop("`x` must not contain Inf or -Inf.", call. = FALSE)
  }
  if (any(is.nan(x))) stop("`x` must not contain NaN.", call. = FALSE)
}

.require_namespace <- function(package, feature) {
  if (!requireNamespace(package, quietly = TRUE)) {
    stop(sprintf("%s requires the optional `%s` package; install it with `install.packages(\"%s\")`.",
                 feature, package, package), call. = FALSE)
  }
}

.validate_qc_report <- function(qc, group_by = NULL) {
  required <- c(".obs", "id", "metric", "value", "lower", "upper", "direction", "is_outlier")
  if (!is.data.frame(qc) || !all(required %in% names(qc))) {
    stop("`qc` must be a `mad_qc()` report containing: .obs, id, metric, value, lower, upper, direction, and is_outlier.", call. = FALSE)
  }
  if (!is.numeric(qc$.obs) || anyNA(qc$.obs) || any(qc$.obs <= 0) || any(qc$.obs != floor(qc$.obs)) ||
      !is.atomic(qc$id) || anyNA(qc$id) ||
      !is.character(qc$metric) || anyNA(qc$metric) || any(qc$metric == "") ||
      !is.numeric(qc$value) || !is.numeric(qc$lower) || !is.numeric(qc$upper) ||
      !is.character(qc$direction) || anyNA(qc$direction) || !all(qc$direction %in% c("lower", "upper", "both")) ||
      !is.logical(qc$is_outlier)) {
    stop("`qc` contains malformed MAD QC columns.", call. = FALSE)
  }
  if (anyDuplicated(qc[c(".obs", "metric")])) {
    stop("`qc` contains duplicate observation-metric records.", call. = FALSE)
  }
  if (nrow(qc) && any(vapply(split(qc$id, qc$.obs), function(x) length(unique(x)) != 1L, logical(1)))) {
    stop("`qc` maps an observation index to multiple identifiers.", call. = FALSE)
  }
  if (!is.null(group_by) && (!is.character(group_by) || !length(group_by) ||
      anyNA(group_by) || any(group_by == "") || anyDuplicated(group_by))) {
    stop("`group_by` must be NULL or one or more unique column names.", call. = FALSE)
  }
  missing_groups <- setdiff(group_by, names(qc))
  if (length(missing_groups)) {
    stop(sprintf("Missing grouping column(s): %s.", paste(missing_groups, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

.mad_stats <- function(x, constant, na_rm, zero_mad) {
  .numeric_vector(x)
  .arg_positive(constant, "constant")
  .arg_flag(na_rm, "na_rm")
  zero_mad <- match.arg(zero_mad, c("zero", "na", "error"))
  if (!na_rm && anyNA(x)) return(list(median = NA_real_, mad = NA_real_, zero = FALSE))
  values <- x[!is.na(x)]
  if (!length(values)) return(list(median = NA_real_, mad = NA_real_, zero = FALSE))
  centre <- stats::median(values)
  spread <- stats::mad(values, center = centre, constant = constant)
  if (spread == 0 && zero_mad == "error") {
    stop("MAD is zero; choose `zero_mad = \"zero\"` or `\"na\"`.", call. = FALSE)
  }
  list(median = centre, mad = spread, zero = spread == 0, zero_mad = zero_mad)
}

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

.numeric_matrix <- function(x) {
  if (inherits(x, "sparseMatrix")) {
    stop("Sparse matrices are not supported; convert explicitly only if safe.", call. = FALSE)
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`x` must be an ordinary numeric matrix.", call. = FALSE)
  }
  if (any(is.infinite(x), na.rm = TRUE)) {
    stop("`x` must not contain Inf or -Inf.", call. = FALSE)
  }
  if (any(is.nan(x))) stop("`x` must not contain NaN.", call. = FALSE)
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

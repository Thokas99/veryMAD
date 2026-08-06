#' Robust MAD scaling
#'
#' Scale a numeric vector, matrix, or numeric data frame by its median and
#' median absolute deviation. For matrices, `margin = 1` scales rows and
#' `margin = 2` scales columns.
#'
#' @param x A numeric vector, matrix, or numeric data frame.
#' @param center Subtract the median?
#' @param scale Divide by the MAD?
#' @param constant Positive MAD consistency constant.
#' @param na_rm Remove missing values when calculating medians and MADs?
#' @param zero_mad Use zero, return `NA`, or error for a zero MAD.
#' @param margin Matrix margin, `1` for rows or `2` for columns. Names
#'   `"rows"` and `"columns"` are also accepted.
#' @return A numeric vector or matrix. Matrix input preserves dimensions and
#'   dimnames; data-frame input returns a matrix.
#' @export
#' @rdname mad_scale
#' @examples
#' mad_scale(c(a = 1, b = 2, c = 100))
#' mad_scale(matrix(1:12, nrow = 3), margin = 1)
mad_scale <- function(x, center = TRUE, scale = TRUE, constant = 1.4826,
                      na_rm = TRUE, zero_mad = c("zero", "na", "error"), margin = 2) {
  .validate_scalar_flag(center, "center"); .validate_scalar_flag(scale, "scale"); .validate_scalar_flag(na_rm, "na_rm")
  .validate_positive(constant, "constant"); zero_mad <- match.arg(zero_mad)
  if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1)))) stop("All data-frame columns must be numeric.", call. = FALSE)
    x <- as.matrix(x)
  }
  if (is.null(dim(x))) return(.scale_vector(x, center, scale, constant, na_rm, zero_mad))
  if (!is.matrix(x) || !is.numeric(x)) stop("`x` must be a numeric matrix or vector.", call. = FALSE)
  if (any(is.infinite(x), na.rm = TRUE)) stop("`x` must not contain Inf or -Inf.", call. = FALSE)
  if (any(is.nan(x))) stop("`x` must not contain NaN.", call. = FALSE)
  if (!length(x)) return(x)
  margin <- if (is.character(margin)) match.arg(margin, c("rows", "columns")) else {
    if (!is.numeric(margin) || length(margin) != 1L || !margin %in% c(1, 2)) stop("`margin` must be 1, 2, `rows`, or `columns`.", call. = FALSE)
    as.integer(margin)
  }
  out <- if (identical(margin, 1L) || identical(margin, "rows")) {
    centres <- matrixStats::rowMedians(x, na.rm = na_rm, useNames = FALSE)
    spreads <- matrixStats::rowMads(x, na.rm = na_rm, constant = constant, useNames = FALSE)
    sweep(.scale_margins(x, centres, spreads, center, scale, zero_mad, 1L), 1L, 0, "+")
  } else {
    centres <- matrixStats::colMedians(x, na.rm = na_rm, useNames = FALSE)
    spreads <- matrixStats::colMads(x, na.rm = na_rm, constant = constant, useNames = FALSE)
    .scale_margins(x, centres, spreads, center, scale, zero_mad, 2L)
  }
  dimnames(out) <- dimnames(x); out
}

.scale_vector <- function(x, center, scale, constant, na_rm, zero_mad) {
  if (!is.numeric(x) || !is.null(dim(x))) stop("`x` must be a numeric vector.", call. = FALSE)
  if (any(is.infinite(x), na.rm = TRUE) || any(is.nan(x))) stop("`x` must not contain Inf, -Inf, or NaN.", call. = FALSE)
  if (!length(x)) return(stats::setNames(numeric(), names(x)))
  centre <- if (na_rm) stats::median(x, na.rm = TRUE) else stats::median(x)
  spread <- if (na_rm) stats::mad(x, center = centre, constant = constant, na.rm = TRUE) else stats::mad(x, center = centre, constant = constant, na.rm = FALSE)
  if (scale && isTRUE(spread == 0)) {
    if (zero_mad == "error") stop("MAD is zero.", call. = FALSE)
    if (zero_mad == "na") return(rep(NA_real_, length(x)))
    spread <- 1
  }
  out <- if (center) x - centre else x
  if (scale) out <- out / spread
  out[is.na(x)] <- NA_real_; names(out) <- names(x); out
}

.scale_margins <- function(x, centres, spreads, center, scale, zero_mad, margin) {
  affected <- which(scale & !is.na(spreads) & spreads == 0)
  if (length(affected) && zero_mad == "error") stop(sprintf("MAD is zero for %s %s.", if (margin == 1) "row" else "column", affected[1]), call. = FALSE)
  safe <- spreads; safe[affected] <- 1
  out <- x
  if (center) out <- sweep(out, margin, centres, "-")
  if (scale) out <- sweep(out, margin, safe, "/")
  if (length(affected) && zero_mad == "na") {
    if (margin == 1) out[affected, ] <- NA_real_ else out[, affected] <- NA_real_
  }
  out[is.na(x)] <- NA_real_; out
}

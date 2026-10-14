#' Center and scale with the median and MAD
#'
#' Optionally center on the median and scale by the median absolute deviation
#' (MAD). For matrices, `margin = 1` operates on rows and `margin = 2` on
#' columns.
#'
#' @param x Numeric vector, matrix, or data frame with numeric columns.
#' @param center Subtract the median from each vector or matrix margin?
#' @param scale Divide each vector or matrix margin by its MAD?
#' @param constant Positive consistency constant applied to the MAD. The
#'   default, 1.4826, makes it consistent with the standard deviation under a
#'   normal distribution.
#' @param na_rm Omit missing values when estimating medians and MADs. Missing
#'   input values remain missing in the result.
#' @param zero_mad Policy for a zero MAD when `scale = TRUE`: `"na"` returns
#'   `NA` for the affected vector or margin, `"zero"` uses a unit denominator,
#'   and `"error"` stops. Centered results under `"zero"` are deviations from
#'   the median without MAD scaling.
#' @param margin Matrix margin: `1` or `"rows"` operates on rows; `2` or
#'   `"columns"` operates on columns. Ignored for vectors.
#' @return A numeric vector or matrix. Matrices preserve dimensions and
#'   dimnames; data frames are returned as matrices.
#' @export
#' @rdname mad_scale
#' @examples
#' gene_expression <- c(gene_A = 8.2, gene_B = 8.7, gene_C = 18.0)
#' mad_scale(gene_expression)
#'
#' expression_matrix <- matrix(
#'   1:12,
#'   nrow = 3,
#'   dimnames = list(paste0("gene_", 1:3), paste0("sample_", 1:4))
#' )
#' mad_scale(expression_matrix, margin = 1)
mad_scale <- function(x, center = TRUE, scale = TRUE, constant = 1.4826,
                      na_rm = TRUE, zero_mad = c("na", "zero", "error"), margin = 2) {
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
  if (!center && !scale) return(x)
  out <- if (identical(margin, 1L) || identical(margin, "rows")) {
    centres <- if (center) matrixStats::rowMedians(x, na.rm = na_rm, useNames = FALSE) else NULL
    spreads <- if (scale) matrixStats::rowMads(x, center = centres, na.rm = na_rm, constant = constant, useNames = FALSE) else NULL
    .scale_margins(x, centres, spreads, center, scale, zero_mad, 1L)
  } else {
    centres <- if (center) matrixStats::colMedians(x, na.rm = na_rm, useNames = FALSE) else NULL
    spreads <- if (scale) matrixStats::colMads(x, center = centres, na.rm = na_rm, constant = constant, useNames = FALSE) else NULL
    .scale_margins(x, centres, spreads, center, scale, zero_mad, 2L)
  }
  dimnames(out) <- dimnames(x); out
}

.scale_vector <- function(x, center, scale, constant, na_rm, zero_mad) {
  if (!is.numeric(x) || !is.null(dim(x))) stop("`x` must be a numeric vector.", call. = FALSE)
  if (any(is.infinite(x), na.rm = TRUE) || any(is.nan(x))) stop("`x` must not contain Inf, -Inf, or NaN.", call. = FALSE)
  if (!length(x)) return(stats::setNames(numeric(), names(x)))
  centre <- if (center) stats::median(x, na.rm = na_rm) else NULL
  spread <- if (scale) {
    if (center) stats::mad(x, center = centre, constant = constant, na.rm = na_rm) else {
      stats::mad(x, constant = constant, na.rm = na_rm)
    }
  } else NULL
  zero <- scale && isTRUE(spread == 0)
  if (zero) {
    if (zero_mad == "error") stop("MAD is zero.", call. = FALSE)
  }
  out <- if (zero && zero_mad == "na") rep(NA_real_, length(x)) else x
  if (center && !(zero && zero_mad == "na")) out <- out - centre
  if (scale && !(zero && zero_mad == "na")) out <- out / if (zero) 1 else spread
  out[is.na(x)] <- NA_real_; names(out) <- names(x); out
}

.scale_margins <- function(x, centres, spreads, center, scale, zero_mad, margin) {
  affected <- if (scale) which(!is.na(spreads) & spreads == 0) else integer()
  if (length(affected) && zero_mad == "error") stop(sprintf("MAD is zero for %s %s.", if (margin == 1) "row" else "column", affected[1]), call. = FALSE)
  safe <- spreads
  if (length(affected)) safe[affected] <- 1
  out <- x
  if (center) out <- sweep(out, margin, centres, "-")
  if (scale) out <- sweep(out, margin, safe, "/")
  if (length(affected) && zero_mad == "na") {
    if (margin == 1) out[affected, ] <- NA_real_ else out[, affected] <- NA_real_
  }
  out[is.na(x)] <- NA_real_; out
}

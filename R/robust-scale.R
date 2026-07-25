#' Median-center and MAD-scale a matrix
#'
#' @param x An ordinary numeric matrix. Sparse matrices are rejected because
#'   median centering can turn implicit zeros into nonzero values and destroy
#'   sparsity; no implicit dense conversion is performed.
#' @param margin Scale across `"rows"` or `"columns"`.
#' @param center Median-center each margin?
#' @param scale MAD-scale each margin?
#' @param constant Positive MAD consistency constant.
#' @param na_rm Remove missing values for summaries?
#' @param zero_mad Behavior for margins with zero MAD.
#' @return A numeric matrix preserving dimensions and dimnames.
#' @export
#' @examples
#' robust_scale(matrix(1:12, nrow = 3))
robust_scale <- function(x, margin = c("rows", "columns"), center = TRUE, scale = TRUE,
                         constant = 1.4826, na_rm = TRUE,
                         zero_mad = c("zero", "na", "error")) {
  .numeric_matrix(x)
  margin <- match.arg(margin); zero_mad <- match.arg(zero_mad)
  .arg_flag(center, "center"); .arg_flag(scale, "scale"); .arg_flag(na_rm, "na_rm")
  .arg_positive(constant, "constant")
  original_dimnames <- dimnames(x)
  work <- if (margin == "columns") t(x) else x
  centres <- if (center) matrixStats::rowMedians(work, na.rm = na_rm) else rep(0, nrow(work))
  spreads <- if (scale) matrixStats::rowMads(work, na.rm = na_rm, constant = constant) else rep(1, nrow(work))
  zero <- scale & !is.na(spreads) & spreads == 0
  if (any(zero) && zero_mad == "error") stop("MAD is zero for at least one matrix margin.", call. = FALSE)
  spreads[zero] <- 1
  out <- sweep(sweep(work, 1L, centres), 1L, spreads, "/")
  if (any(zero)) out[zero, ] <- if (zero_mad == "zero") 0 else NA_real_
  out[is.na(work)] <- NA_real_
  if (margin == "columns") out <- t(out)
  dimnames(out) <- original_dimnames
  out
}

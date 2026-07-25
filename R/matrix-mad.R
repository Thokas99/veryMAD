#' Row- and column-wise MAD
#'
#' @param x An ordinary numeric matrix. Sparse matrices are rejected.
#' @param na_rm Remove missing values?
#' @param constant Positive MAD consistency constant.
#' @return A named numeric vector, when the summarized margin has names.
#' @export
#' @examples
#' m <- matrix(1:12, nrow = 3)
#' row_mad(m)
row_mad <- function(x, na_rm = TRUE, constant = 1.4826) {
  .numeric_matrix(x); .arg_flag(na_rm, "na_rm"); .arg_positive(constant, "constant")
  out <- matrixStats::rowMads(x, na.rm = na_rm, constant = constant)
  names(out) <- rownames(x); out
}

#' @rdname row_mad
#' @export
col_mad <- function(x, na_rm = TRUE, constant = 1.4826) {
  .numeric_matrix(x); .arg_flag(na_rm, "na_rm"); .arg_positive(constant, "constant")
  out <- matrixStats::colMads(x, na.rm = na_rm, constant = constant)
  names(out) <- colnames(x); out
}

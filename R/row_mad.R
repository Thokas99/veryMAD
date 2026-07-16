#' Calculate MAD across rows or columns of a matrix
#'
#' Computes the median absolute deviation for each row (or column) of a
#' numeric matrix, e.g. to rank variable genes, proteins, or other features
#' in an omics matrix. Uses [matrixStats::rowMads()]/[matrixStats::colMads()]
#' internally for speed.
#'
#' @param x A numeric matrix.
#' @param margin `1` to compute MAD across rows (the default; e.g. one row
#'   per feature), or `2` to compute across columns.
#' @param na_rm Logical. Should missing values be removed?
#' @param constant A positive scaling constant passed to
#'   [matrixStats::rowMads()]/[matrixStats::colMads()].
#' @param verbose Logical. If `TRUE`, report the matrix margin summarized.
#'
#' @return A numeric vector of MAD values, one per row (or column) of `x`,
#'   named with the corresponding row (or column) names if present.
#' @export
#'
#' @family robust MAD helpers
#'
#' @examples
#' m <- matrix(c(1, 2, 2, 3, 100, 10, 20, 20, 30, 40), nrow = 2, byrow = TRUE)
#' row_mad(m)
row_mad <- function(x, margin = 1, na_rm = FALSE, constant = 1.4826, verbose = FALSE) {
  if (!is.matrix(x) || !is.numeric(x)) {
    cli::cli_abort("{.arg x} must be a numeric matrix.")
  }

  if (!margin %in% c(1L, 2L)) {
    cli::cli_abort("{.arg margin} must be 1 (rows) or 2 (columns).")
  }

  if (margin == 1L) {
    result <- matrixStats::rowMads(x, na.rm = na_rm, constant = constant)
    names(result) <- rownames(x)
  } else {
    result <- matrixStats::colMads(x, na.rm = na_rm, constant = constant)
    names(result) <- colnames(x)
  }

  if (isTRUE(verbose)) {
    cli::cli_inform("Calculated MAD for {length(result)} matrix {if (margin == 1L) 'row' else 'column'}{?s}.")
  }
  result
}

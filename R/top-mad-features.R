#' Rank features by raw marginal MAD
#'
#' This is not mean-variance-aware and is not a replacement for Seurat or
#' Scanpy highly variable feature selection. It may perform poorly on raw
#' sparse counts.
#' @param x An ordinary numeric matrix with features in rows.
#' @param n Number of feature names or indices to return.
#' @inheritParams row_mad
#' @return Row names, or row indices, ordered by decreasing raw MAD.
#' @export
#' @examples
#' top_mad_features(matrix(rnorm(50), nrow = 10), n = 3)
top_mad_features <- function(x, n = 2000, na_rm = TRUE, constant = 1.4826) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n <= 0 || n != as.integer(n)) {
    stop("`n` must be one finite positive whole number.", call. = FALSE)
  }
  scores <- row_mad(x, na_rm, constant)
  labels <- rownames(x) %||% seq_along(scores)
  order <- order(scores, decreasing = TRUE, na.last = TRUE)
  labels[order[seq_len(min(n, length(scores)))]]
}

`%||%` <- function(x, y) if (is.null(x)) y else x

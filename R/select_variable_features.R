#' Select the most variable features by MAD
#'
#' Ranks features (rows, or columns if `margin = 2`) of a matrix by
#' variability and returns the top `n`. A thin convenience wrapper around
#' [row_mad()].
#'
#' @param x A numeric matrix.
#' @param n The number of features to select.
#' @param margin `1` to rank rows (the default; e.g. one row per feature),
#'   or `2` to rank columns.
#' @param na_rm Logical. Should missing values be removed?
#' @param constant A positive scaling constant passed to [row_mad()].
#'
#' @return A character vector of the top `n` row (or column) names, ordered
#'   by decreasing MAD. If `x` has no names on that margin, an integer
#'   vector of indices is returned instead.
#' @export
#'
#' @examples
#' m <- matrix(rnorm(50), nrow = 10)
#' rownames(m) <- paste0("gene", 1:10)
#' select_variable_features(m, n = 3)
select_variable_features <- function(
  x,
  n = 2000,
  margin = 1,
  na_rm = FALSE,
  constant = 1.4826
) {
  if (!is.numeric(n) || length(n) != 1L || n <= 0) {
    cli::cli_abort("{.arg n} must be one positive number.")
  }

  scores <- row_mad(x, margin = margin, na_rm = na_rm, constant = constant)

  labels <- names(scores)
  if (is.null(labels)) {
    labels <- seq_along(scores)
  }

  n <- min(n, length(scores))
  labels[order(scores, decreasing = TRUE)][seq_len(n)]
}

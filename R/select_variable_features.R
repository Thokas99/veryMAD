#' Select the most variable features by MAD
#'
#' Ranks features (rows, or columns if `margin = 2`) of a matrix by
#' variability and returns the top `n`. A thin convenience wrapper around
#' [row_mad()].
#'
#' @param x A numeric matrix.
#' @param n A finite positive whole number of features to select.
#' @param margin `1` to rank rows (the default; e.g. one row per feature),
#'   or `2` to rank columns.
#' @param na_rm Logical. Should missing values be removed?
#' @param constant A positive scaling constant passed to [row_mad()].
#' @param verbose Logical. If `TRUE`, report how many features were selected.
#'
#' @return A character vector of the top `n` row (or column) names, ordered
#'   by decreasing MAD. If `x` has no names on that margin, an integer
#'   vector of indices is returned instead.
#' @export
#'
#' @family robust MAD helpers
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
  constant = 1.4826,
  verbose = FALSE
) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n <= 0 || n != as.integer(n)) {
    cli::cli_abort("{.arg n} must be one finite positive whole number.")
  }

  scores <- row_mad(x, margin = margin, na_rm = na_rm, constant = constant)

  labels <- names(scores)
  if (is.null(labels)) {
    labels <- seq_along(scores)
  }

  n <- min(n, length(scores))
  result <- labels[order(scores, decreasing = TRUE)][seq_len(n)]
  if (isTRUE(verbose)) {
    cli::cli_inform("Selected {length(result)} variable feature{?s} by MAD.")
  }
  result
}

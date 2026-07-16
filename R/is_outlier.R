#' Flag robust outliers
#'
#' Flags values whose robust MAD score (see [mad_score()]) exceeds a
#' threshold in absolute value.
#'
#' @param x A numeric vector.
#' @param threshold A positive number. Values with `abs(mad_score(x)) >
#'   threshold` are flagged as outliers. The default of 3.5 follows the
#'   modified z-score convention of Iglewicz and Hoaglin (1993).
#' @param constant A positive scaling constant passed to [stats::mad()].
#' @param na_rm Logical. Should missing values be removed?
#'
#' @return A logical vector the same length as `x`. `NA` where `x` is
#'   missing.
#' @export
#'
#' @examples
#' is_outlier(c(1, 2, 2, 3, 100))
is_outlier <- function(x, threshold = 3.5, constant = 1.4826, na_rm = FALSE) {
  if (!is.numeric(threshold) || length(threshold) != 1L || threshold <= 0) {
    cli::cli_abort("{.arg threshold} must be one positive number.")
  }

  scores <- mad_score(x, constant = constant, na_rm = na_rm)
  abs(scores) > threshold
}

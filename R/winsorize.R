#' Winsorize outliers using a robust MAD threshold
#'
#' Caps values whose robust MAD score (see [mad_score()]) exceeds a
#' threshold in absolute value, replacing them with the nearest boundary
#' value rather than removing them.
#'
#' @param x A numeric vector.
#' @param threshold A positive number. Values with `abs(mad_score(x)) >
#'   threshold` are capped at the corresponding boundary.
#' @param constant A positive scaling constant passed to [stats::mad()].
#' @param na_rm Logical. Should missing values be removed?
#'
#' @return A numeric vector the same length as `x`, with outliers capped.
#' @export
#'
#' @examples
#' winsorize(c(1, 2, 2, 3, 100))
winsorize <- function(x, threshold = 3.5, constant = 1.4826, na_rm = FALSE) {
  if (!is.numeric(threshold) || length(threshold) != 1L || threshold <= 0) {
    cli::cli_abort("{.arg threshold} must be one positive number.")
  }

  centre <- stats::median(x, na.rm = na_rm)
  scale <- stats::mad(x, center = centre, constant = constant, na.rm = na_rm)

  if (is.na(scale) || scale == 0) {
    return(x)
  }

  lower <- centre - threshold * scale
  upper <- centre + threshold * scale

  pmin(pmax(x, lower), upper)
}

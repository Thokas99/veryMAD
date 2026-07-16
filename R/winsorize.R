#' Winsorize outliers using a robust MAD threshold
#'
#' Caps values whose robust MAD score (see [mad_score()]) exceeds a
#' threshold in absolute value, replacing them with the nearest boundary
#' value rather than removing them.
#'
#' @param x A numeric vector.
#' @param threshold A finite positive number. Values with `abs(mad_score(x)) >
#'   threshold` are capped at the corresponding boundary.
#' @param constant A positive scaling constant passed to [stats::mad()].
#' @param na_rm Logical. Should missing values be removed?
#' @param verbose Logical. If `TRUE`, report how many values were capped.
#'
#' @return A numeric vector the same length as `x`, with outliers capped.
#' @export
#'
#' @family robust MAD helpers
#'
#' @examples
#' winsorize(c(1, 2, 2, 3, 100))
winsorize <- function(x, threshold = 3.5, constant = 1.4826, na_rm = FALSE, verbose = FALSE) {
  if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold) || threshold <= 0) {
    cli::cli_abort("{.arg threshold} must be one finite positive number.")
  }

  centre <- stats::median(x, na.rm = na_rm)
  scale <- stats::mad(x, center = centre, constant = constant, na.rm = na_rm)

  if (is.na(scale) || scale == 0) {
    return(x)
  }

  lower <- centre - threshold * scale
  upper <- centre + threshold * scale

  result <- pmin(pmax(x, lower), upper)
  if (isTRUE(verbose)) {
    cli::cli_inform("Winsorized {sum(result != x, na.rm = TRUE)} value{?s} using MAD threshold {threshold}.")
  }
  result
}

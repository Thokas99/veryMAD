#' Calculate robust MAD scores
#'
#' Calculates standardized deviations from the median using the median
#' absolute deviation.
#'
#' @param x A numeric vector.
#' @param constant A positive scaling constant passed to [stats::mad()].
#' @param na_rm Logical. Should missing values be removed?
#' @param verbose Logical. If `TRUE`, report a short summary with [cli::cli_inform()].
#'
#' @return A numeric vector with one robust score per input value.
#' @export
#'
#' @family robust MAD helpers
#'
#' @examples
#' mad_score(c(1, 2, 2, 3, 100))
mad_score <- function(x, constant = 1.4826, na_rm = FALSE, verbose = FALSE) {
  if (!is.numeric(x)) {
    cli::cli_abort("{.arg x} must be a numeric vector.")
  }

  if (!is.numeric(constant) || length(constant) != 1L || constant <= 0) {
    cli::cli_abort("{.arg constant} must be one positive number.")
  }

  centre <- stats::median(x, na.rm = na_rm)
  scale <- stats::mad(
    x,
    center = centre,
    constant = constant,
    na.rm = na_rm
  )

  if (is.na(scale) || scale == 0) {
    return(rep(NA_real_, length(x)))
  }

  result <- (x - centre) / scale
  if (isTRUE(verbose)) {
    cli::cli_inform("Calculated MAD scores for {length(x)} value{?s}.")
  }
  result
}

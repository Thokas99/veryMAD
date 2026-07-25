#' Calculate MAD-based limits
#'
#' @param x A numeric vector without infinite values.
#' @param nmads Number of MADs from the median.
#' @param direction Tail(s) for which to calculate finite limits.
#' @param constant Positive MAD consistency constant.
#' @param na_rm Remove missing values before calculation?
#' @param zero_mad Behavior when MAD is zero: neutral limits, missing limits,
#'   or an error.
#' @return A one-row tibble with `median`, `mad`, `lower`, and `upper`.
#' @export
#' @examples
#' mad_limits(c(1, 2, 2, 3, 100), direction = "upper")
mad_limits <- function(x, nmads = 3, direction = c("both", "lower", "upper"),
                       constant = 1.4826, na_rm = TRUE,
                       zero_mad = c("zero", "na", "error")) {
  direction <- match.arg(direction)
  zero_mad <- match.arg(zero_mad)
  .arg_positive(nmads, "nmads")
  s <- .mad_stats(x, constant, na_rm, zero_mad)
  limits <- .mad_limits_from_stats(s, nmads, direction, zero_mad)
  tibble::tibble(median = s$median, mad = s$mad,
    lower = unname(limits[["lower"]]), upper = unname(limits[["upper"]]))
}

#' Calculate MAD-based limits
#'
#' @param x A numeric vector without infinite values.
#' @param nmads Number of MADs from the median.
#' @param direction Tail(s) for which to calculate finite limits.
#' @param constant Positive MAD consistency constant.
#' @param na_rm Remove missing values before calculation?
#' @param zero_mad Behavior when MAD is zero: neutral limits, missing limits,
#'   or an error.
#' @return A one-row data frame with `median`, `mad`, `lower`, and `upper`.
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
  lower <- if (direction %in% c("both", "lower")) s$median - nmads * s$mad else -Inf
  upper <- if (direction %in% c("both", "upper")) s$median + nmads * s$mad else Inf
  if (isTRUE(s$zero)) {
    if (zero_mad == "zero") { lower <- -Inf; upper <- Inf }
    if (zero_mad == "na") { lower <- NA_real_; upper <- NA_real_ }
  }
  data.frame(median = s$median, mad = s$mad, lower = lower, upper = upper)
}

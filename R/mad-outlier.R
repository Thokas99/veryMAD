#' Flag MAD outliers
#'
#' @inheritParams mad_limits
#' @return A logical vector with the length and names of `x`.
#' @export
#' @examples
#' is_mad_outlier(c(1, 2, 2, 3, 100), direction = "upper")
is_mad_outlier <- function(x, nmads = 3, direction = c("both", "lower", "upper"),
                           constant = 1.4826, na_rm = TRUE,
                           zero_mad = c("zero", "na", "error")) {
  direction <- match.arg(direction)
  zero_mad <- match.arg(zero_mad)
  .arg_positive(nmads, "nmads")
  stats <- .mad_stats(x, constant, na_rm, zero_mad)
  .mad_flags(x, stats, nmads, direction, zero_mad)
}

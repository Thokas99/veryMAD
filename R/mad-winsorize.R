#' Winsorize values at MAD limits
#'
#' @inheritParams mad_limits
#' @return A numeric vector with the length and names of `x`.
#' @export
#' @examples
#' winsorize_mad(c(1, 2, 2, 3, 100), direction = "upper")
winsorize_mad <- function(x, nmads = 3, direction = c("both", "lower", "upper"),
                          constant = 1.4826, na_rm = TRUE,
                          zero_mad = c("zero", "na", "error")) {
  limits <- mad_limits(x, nmads, match.arg(direction), constant, na_rm, match.arg(zero_mad))
  pmin(pmax(x, limits$lower), limits$upper)
}

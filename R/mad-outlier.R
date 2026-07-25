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
  limits <- mad_limits(x, nmads, match.arg(direction), constant, na_rm, match.arg(zero_mad))
  x < limits$lower | x > limits$upper
}

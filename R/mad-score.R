#' Calculate robust MAD scores
#'
#' @inheritParams mad_limits
#' @return A numeric vector with the length and names of `x`.
#' @export
#' @examples
#' mad_score(c(a = 1, b = 2, c = 100))
mad_score <- function(x, constant = 1.4826, na_rm = TRUE,
                      zero_mad = c("zero", "na", "error")) {
  zero_mad <- match.arg(zero_mad)
  s <- .mad_stats(x, constant, na_rm, zero_mad)
  out <- (x - s$median) / s$mad
  if (isTRUE(s$zero)) out[!is.na(x)] <- if (zero_mad == "zero") 0 else NA_real_
  out
}

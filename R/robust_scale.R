#' Median-center and MAD-scale a matrix
#'
#' Robustly standardizes each row (or column) of a matrix by subtracting the
#' median and dividing by the median absolute deviation, analogous to
#' [base::scale()] but resistant to outliers. Well suited to expression
#' matrices, where rows are typically features (genes, proteins, miRNAs) and
#' columns are samples.
#'
#' @param x A numeric matrix.
#' @param margin `1` to scale across rows (the default; e.g. one row per
#'   feature), or `2` to scale across columns.
#' @param center Logical. Should each row/column be median-centered?
#' @param scale Logical. Should each row/column be MAD-scaled?
#' @param zero_mad What to do with rows/columns whose MAD is 0 (common in
#'   sparse omics data, where a feature can be constant within a group).
#'   One of `"zero"` (set scaled values to 0, the default), `"na"` (set
#'   scaled values to `NA`), or `"error"` (abort).
#'
#' @return A numeric matrix the same shape as `x`.
#' @export
#'
#' @family robust MAD helpers
#'
#' @examples
#' m <- matrix(c(1, 2, 2, 3, 100, 10, 20, 20, 30, 40), nrow = 2, byrow = TRUE)
#' robust_scale(m)
robust_scale <- function(
  x,
  margin = 1,
  center = TRUE,
  scale = TRUE,
  zero_mad = c("zero", "na", "error")
) {
  zero_mad <- match.arg(zero_mad)

  if (!is.matrix(x) || !is.numeric(x)) {
    cli::cli_abort("{.arg x} must be a numeric matrix.")
  }

  if (!margin %in% c(1L, 2L)) {
    cli::cli_abort("{.arg margin} must be 1 (rows) or 2 (columns).")
  }

  if (margin == 2L) {
    x <- t(x)
  }

  centres <- if (center) {
    matrixStats::rowMedians(x, na.rm = TRUE)
  } else {
    rep(0, nrow(x))
  }
  mads <- if (scale) {
    matrixStats::rowMads(x, na.rm = TRUE)
  } else {
    rep(1, nrow(x))
  }

  zero_rows <- scale & mads == 0
  if (any(zero_rows)) {
    if (zero_mad == "error") {
      cli::cli_abort(
        "{sum(zero_rows)} row(s) have MAD = 0; cannot scale."
      )
    }
    # Avoid division by zero; overwritten below per `zero_mad`.
    mads[zero_rows] <- 1
  }

  scaled <- (x - centres) / mads

  if (any(zero_rows)) {
    scaled[zero_rows, ] <- switch(
      zero_mad,
      zero = 0,
      na = NA_real_
    )
  }

  if (margin == 2L) {
    scaled <- t(scaled)
  }

  scaled
}

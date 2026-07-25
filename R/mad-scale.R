#' Robust MAD scaling for vectors and dense matrices
#'
#' `mad_scale()` median-centres and MAD-scales a numeric vector or each margin
#' of an ordinary dense numeric matrix. For expression heatmaps, use
#' `margin = "rows"` and pass the result to the plotting function without
#' applying another scaling step.
#'
#' @param x A numeric vector or ordinary dense numeric matrix.
#' @param margin Matrix margin to scale. Ignored for vector input.
#' @param constant Positive MAD consistency constant.
#' @param na_rm Remove missing values when calculating centres and spreads?
#' @param zero_mad Behavior for zero-MAD vectors or matrix margins.
#' @return A numeric vector or matrix preserving names or dimnames.
#' @export
#' @examples
#' mad_scale(c(a = 1, b = 2, c = 100))
#' x <- matrix(1:12, nrow = 3, dimnames = list(paste0("gene", 1:3), NULL))
#' mad_scale(x, margin = "rows")
mad_scale <- function(x, margin = c("rows", "columns"), constant = 1.4826,
                      na_rm = TRUE, zero_mad = c("zero", "na", "error")) {
  if (is.null(dim(x))) {
    return(mad_score(x, constant = constant, na_rm = na_rm, zero_mad = zero_mad))
  }
  .numeric_matrix(x)
  margin <- match.arg(margin)
  zero_mad <- match.arg(zero_mad)
  .arg_positive(constant, "constant")
  .arg_flag(na_rm, "na_rm")
  storage.mode(x) <- "double"
  if (margin == "rows") {
    centers <- matrixStats::rowMedians(x, na.rm = na_rm, useNames = FALSE)
    spreads <- matrixStats::rowMads(x, na.rm = na_rm, constant = constant, useNames = FALSE)
    out <- sweep(x, 1L, centers, FUN = "-")
    affected <- which(!is.na(spreads) & spreads == 0)
    safe_spreads <- spreads
    safe_spreads[affected] <- 1
    out <- sweep(out, 1L, safe_spreads, FUN = "/")
    out <- .handle_zero_margins(out, x, affected, 1L, zero_mad)
  } else {
    centers <- matrixStats::colMedians(x, na.rm = na_rm, useNames = FALSE)
    spreads <- matrixStats::colMads(x, na.rm = na_rm, constant = constant, useNames = FALSE)
    out <- sweep(x, 2L, centers, FUN = "-")
    affected <- which(!is.na(spreads) & spreads == 0)
    safe_spreads <- spreads
    safe_spreads[affected] <- 1
    out <- sweep(out, 2L, safe_spreads, FUN = "/")
    out <- .handle_zero_margins(out, x, affected, 2L, zero_mad)
  }
  out[is.na(x)] <- NA_real_
  dimnames(out) <- dimnames(x)
  out
}

.handle_zero_margins <- function(out, x, affected, margin, zero_mad) {
  if (!length(affected)) return(out)
  if (zero_mad == "error") {
    labels <- if (margin == 1L) rownames(x) else colnames(x)
    label <- if (is.null(labels) || is.na(labels[affected[1L]]) || labels[affected[1L]] == "") {
      as.character(affected[1L])
    } else {
      sprintf("`%s`", labels[affected[1L]])
    }
    stop(sprintf("MAD is zero for %s %s.", if (margin == 1L) "row" else "column", label), call. = FALSE)
  }
  value <- if (zero_mad == "zero") 0 else NA_real_
  if (margin == 1L) out[affected, ] <- value else out[, affected] <- value
  out
}

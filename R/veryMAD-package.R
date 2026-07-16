#' veryMAD: Robust MAD-based QC helpers
#'
#' veryMAD provides small helpers for robust outlier detection, scaling, and
#' variable feature selection using the median absolute deviation (MAD).
#'
#' The package is designed for quality-control metrics and expression-like
#' matrices where a few extreme values should not drive the center or scale.
#'
#' @section Main functions:
#' - [mad_score()] returns robust z-like scores based on the median and MAD.
#' - [is_outlier()] flags observations whose absolute MAD score exceeds a
#'   threshold.
#' - [winsorize()] caps outliers at MAD-based limits.
#' - [row_mad()] computes row-wise or column-wise MAD values for numeric
#'   matrices.
#' - [robust_scale()] centers and scales numeric matrices with medians and MADs.
#' - [mad_qc()] adds MAD-based QC outlier flags to tabular QC metrics.
#' - [select_variable_features()] selects matrix rows or columns with the
#'   largest MAD values.
#'
#' @keywords internal
"_PACKAGE"

#' Run explicit MAD QC on Seurat cell metadata
#'
#' `mad_qc_seurat()` preserves the v0.3.2 long report and per-metric annotation
#' contract. New code should use [mad_qc_sc()] for compact reports and one-column
#' annotation. Neither function filters cells.
#'
#' Use explicit transformations for count-like cell metrics when needed. For
#' example, `nCount_RNA` and `nFeature_RNA` are commonly right-skewed and can be
#' thresholded on `log1p`, while bounded percentages such as `percent.mt` usually
#' remain on the raw scale. Upper-tail count or feature flags are inspection
#' warnings, not doublet calls.
#'
#' @param object A Seurat object.
#' @param metrics A named character vector mapping metadata columns to
#'   `"lower"`, `"upper"`, or `"both"`.
#' @param nmads Number of MADs from the median used to define thresholds.
#' @param transform Optional named character vector of explicit per-metric
#'   transformations passed to [mad_qc()].
#' @param group_by Optional metadata column used to estimate thresholds within
#'   groups.
#' @param action Return a QC `"report"` or `"annotate"` the object metadata.
#' @param overwrite Allow annotation columns to replace existing columns.
#' @param verbose Print the concise QC summary?
#' @param zero_mad Compatibility policy for groups with zero MAD: `"zero"`,
#'   `"na"`, or `"error"`.
#'
#' @return For `action = "report"`, a `mad_qc` data frame. For
#'   `action = "annotate"`, the input Seurat object with metadata flag columns
#'   added.
#' @export
#'
#' @examples
#' metrics <- c(
#'   nCount_RNA = "both",
#'   nFeature_RNA = "both",
#'   percent.mt = "upper"
#' )
#' transform <- c(
#'   nCount_RNA = "log1p",
#'   nFeature_RNA = "log1p"
#' )
#' # report <- mad_qc_seurat(seurat_object, metrics, transform = transform)
#' # seurat_object <- mad_qc_seurat(
#' #   seurat_object, metrics, transform = transform, action = "annotate"
#' # )
mad_qc_seurat <- function(object, metrics = NULL, nmads = 3, transform = NULL,
                          group_by = NULL, zero_mad = c("zero", "na", "error"),
                          action = c("report", "annotate"), overwrite = FALSE,
                          verbose = FALSE) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    cli::cli_abort(c("x" = "{.fn mad_qc_seurat} requires the optional {.pkg SeuratObject} package."))
  }
  if (!inherits(object, "Seurat")) cli::cli_abort(c("x" = "{.arg object} must be a Seurat object."))
  metrics <- .validate_wrapper_metrics(object[[]], metrics, "mad_qc_seurat")
  action <- match.arg(action); zero_mad <- match.arg(zero_mad); .arg_flag(overwrite, "overwrite"); .arg_flag(verbose, "verbose")
  qc <- mad_qc(object[[]], metrics, nmads, group_by, transform, zero_mad = zero_mad)
  if (verbose) {
    compact <- .mad_qc_compact(object[[]], metrics, nmads,
      if (is.null(transform)) "none" else transform, group_by)
    .inform_mad_qc(compact, "single-cell QC", "cells")
  }
  if (action == "report") return(qc)
  annotated <- annotate_mad_qc(object[[]], qc, overwrite = overwrite)
  additions <- annotated[c(paste0(names(metrics), "_mad_outlier"), "mad_qc_outlier")]
  rownames(additions) <- rownames(object[[]])
  SeuratObject::AddMetaData(object, metadata = additions)
}

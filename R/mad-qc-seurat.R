#' Run explicit MAD QC on Seurat cell metadata
#'
#' `mad_qc_seurat()` applies [mad_qc()] to `object[[]]` and either returns the
#' long QC report or annotates Seurat metadata with outlier flags. It does not
#' filter cells and does not alter expression data.
#'
#' Use explicit transformations for count-like cell metrics when needed. For
#' example, `nCount_RNA` and `nFeature_RNA` are commonly right-skewed and can be
#' thresholded on `log10p`, while bounded percentages such as `percent.mt` usually
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
#' @param zero_mad How to handle groups with zero MAD.
#' @param action Return a QC `"report"` or `"annotate"` the object metadata.
#' @param overwrite Allow annotation columns to replace existing columns.
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
#'   nCount_RNA = "log10p",
#'   nFeature_RNA = "log10p"
#' )
#' # report <- mad_qc_seurat(seurat_object, metrics, transform = transform)
#' # seurat_object <- mad_qc_seurat(
#' #   seurat_object, metrics, transform = transform, action = "annotate"
#' # )
mad_qc_seurat <- function(object,
                          metrics = c(nCount_RNA = "lower", nFeature_RNA = "lower", percent.mt = "upper"),
                          nmads = 3, transform = NULL, group_by = NULL,
                          zero_mad = c("zero", "na", "error"),
                          action = c("report", "annotate"), overwrite = FALSE) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("`mad_qc_seurat()` requires SeuratObject; install it with `install.packages(\"SeuratObject\")`.", call. = FALSE)
  }
  if (!inherits(object, "Seurat")) stop("`object` must be a Seurat object.", call. = FALSE)
  action <- match.arg(action); zero_mad <- match.arg(zero_mad); .arg_flag(overwrite, "overwrite")
  metadata <- object[[]]
  qc <- mad_qc(metadata, metrics, nmads, group_by, transform, zero_mad = zero_mad)
  if (action == "report") return(qc)
  annotated <- annotate_mad_qc(metadata, qc, overwrite = overwrite)
  additions <- annotated[c(paste0(names(metrics), "_mad_outlier"), "mad_qc_outlier")]
  rownames(additions) <- rownames(metadata)
  SeuratObject::AddMetaData(object, metadata = additions)
}

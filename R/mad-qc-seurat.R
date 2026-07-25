#' Apply MAD QC to Seurat cell metadata
#'
#' @param object A `Seurat` object.
#' @inheritParams mad_qc
#' @param metrics Named QC metrics and directions.
#' @param action Return an annotated object or tidy report.
#' @return A tidy report or annotated Seurat object.
#' @export
mad_qc_seurat <- function(object,
                          metrics = c(nCount_RNA = "lower", nFeature_RNA = "lower", percent.mt = "upper"),
                          nmads = 3, transform = NULL, group_by = NULL,
                          zero_mad = c("zero", "na", "error"),
                          action = c("annotate", "report")) {
  if (!requireNamespace("SeuratObject", quietly = TRUE)) {
    stop("`mad_qc_seurat()` requires SeuratObject; install it with `install.packages(\"SeuratObject\")`.", call. = FALSE)
  }
  if (!inherits(object, "Seurat")) stop("`object` must be a Seurat object.", call. = FALSE)
  action <- match.arg(action); zero_mad <- match.arg(zero_mad)
  metadata <- object[[]]
  qc <- mad_qc(metadata, metrics, nmads, group_by, transform, zero_mad = zero_mad)
  if (action == "report") return(qc)
  annotated <- annotate_mad_qc(metadata, qc)
  additions <- annotated[setdiff(names(annotated), names(metadata))]
  rownames(additions) <- rownames(metadata)
  SeuratObject::AddMetaData(object, metadata = additions)
}

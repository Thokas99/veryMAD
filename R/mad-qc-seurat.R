#' Apply MAD QC to Seurat cell metadata
#'
#' @param object A `Seurat` object.
#' @inheritParams mad_qc
#' @param metrics Named QC metrics and directions.
#' @param action Return a tidy report or annotated object.
#' @param overwrite Replace existing veryMAD flag columns when annotating?
#' @return A tidy report or annotated Seurat object.
#' @export
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

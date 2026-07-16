.onLoad <- function(libname, pkgname) {
  register_seurat_method()
  setHook(
    packageEvent("Seurat", "onLoad"),
    function(...) register_seurat_method()
  )
}

register_seurat_method <- function() {
  if (requireNamespace("Seurat", quietly = TRUE)) {
    S7::method(mad_qc, methods::getClass("Seurat")) <- mad_qc_seurat_method
  }
}

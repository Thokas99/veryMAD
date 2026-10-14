utils::globalVariables(c("label", "metric", "status", "value"))

.verymad_banner <- function(version = as.character(utils::packageVersion("veryMAD"))) {
  paste0("veryMAD ", version, ": MAD quality control and robust scaling")
}

.verymad_startup_message <- function(quiet = getOption("veryMAD.quiet", FALSE),
                                     version = as.character(utils::packageVersion("veryMAD"))) {
  if (!isTRUE(quiet)) packageStartupMessage(.verymad_banner(version))
  invisible(NULL)
}

.onAttach <- function(libname, pkgname) {
  if (interactive()) .verymad_startup_message()
}

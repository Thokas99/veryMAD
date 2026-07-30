.verymad_banner <- function(version = as.character(utils::packageVersion("veryMAD"))) {
  paste(
    "                       __  __    _    ____",
    "__   _____ _ __ _   _|  \\/  |  / \\  |  _ \\",
    "\\ \\ / / _ \\ '__| | | | |\\/| | / _ \\ | | | |",
    " \\ V /  __/ |  | |_| | |  | |/ ___ \\| |_| |",
    "  \\_/ \\___|_|   \\__, |_|  |_/_/   \\_\\____/",
    "                |___/",
    "",
    paste0("veryMAD ", version, " \u2014 robust QC, explicitly yours."),
    sep = "\n"
  )
}

.verymad_startup_message <- function(quiet = getOption("veryMAD.quiet", FALSE),
                                     version = as.character(utils::packageVersion("veryMAD"))) {
  if (!isTRUE(quiet)) packageStartupMessage(.verymad_banner(version))
  invisible(NULL)
}

.onAttach <- function(libname, pkgname) {
  if (interactive()) .verymad_startup_message()
}

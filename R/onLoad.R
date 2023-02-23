.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "Note: Some functions in this package require pre-computed tables of data",
    " which are available externally at github.com/uni-Arya/fdpbandsdata.",
    " The size of the data is roughly 81Mb. You may run",
    " devtools::install_github(\"uni-Arya/fdpbandsdata\") to download these",
    " tables. Make sure to restart R/RStudio after installation.",
    " You can always remove.packages(\"fdpbandsdata\") to delete the data",
    " after installation."
  )
}

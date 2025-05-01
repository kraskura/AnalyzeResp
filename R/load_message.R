.onAttach <- function(libname, pkgname) {
  packageStartupMessage(pkgname, " was last updated: ",  Sys.time())
}

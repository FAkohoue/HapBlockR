#' @useDynLib HapBlockR, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom digest digest
NULL

.onUnload <- function(libpath) {
  library.dynam.unload("HapBlockR", libpath)
}

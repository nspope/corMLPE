#' @import RcppArmadillo
#' @import nlme
#' @import MASS
#' @useDynLib corMLPE
NULL

.onUnload <- function (libpath) {
  library.dynam.unload("corMLPE", libpath)
}


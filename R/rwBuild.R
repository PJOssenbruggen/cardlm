#' \code{rwBuild} estimate parameters of \code{rw} dlm model.
#'
#' @param x, a vector.
#' @usage rwBuild(x)
#' @export
rwBuild <- function(x) {
  return(dlm::dlmModPoly(order = 1, dV = exp(x[1]), dW = exp(x[2])))
}

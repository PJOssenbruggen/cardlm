#' \code{xuBuild} estimate parameters of location \code{x} and speed \code{u} dlm model.
#'
#' @param x, a vector.
#' @usage xuBuild(x)
#' @export
xuBuild <- function(x) {
  Vsd    <- exp(x[1:2])
  Vcorr  <- tanh(x[3])
  V      <- Vsd %o% Vsd
  V[1,2] <- V[2,1] <- V[1,2] * Vcorr
  Wsd    <- exp(x[4:5])
  Wcorr  <- tanh(x[6])
  W      <- Wsd %o% Wsd
  W[1,2] <- W[2,1] <- W[1,2] * Wcorr
  return(list(
    m0 = rep(0,2),
    C0 = 1e7 * diag(2),
    FF = diag(2),
    GG = diag(2),
    V = V,
    W = W))
}

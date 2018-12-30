#' \code{linAcc} is a wrapper function for estimating a successful, safe merge.
#'
#' @usage uses \code{xabparam}, \code{xab} and \code{uab} from the \code{cartools} Package.
# #' @examples
# #' linAcc()
#' @export
linAcc <- function() {
  y        <- ts(rnorm(30,60,10))
  buildMod <- function(x) {
    F.     <- matrix(c(1,0.5,0.125), nrow = 1)
    V      <- matrix(rep(0,9), ncol = 3)
    V[1,1] <- exp(x[1])
    V[2,2] <- exp(x[2])
    V[3,3] <- exp(x[3])
    V[1,2] <- exp(x[4])
    V[2,3] <- exp(x[5])
    m0     <- matrix(60)
    C0     <- matrix(1000)
    G      <- matrix(1)
    W      <- matrix(exp(x[6]))
    Mod    <- list(FF = F., V  = V, GG = G, W  = W, m0 = m0, C0 = C0)
    return(Mod)
  }
  fitMod <- dlmMLE(y, parm = rep(0,6), build = buildMod, debug = TRUE,
                   hessian = TRUE, control = list(maxit = 500))
  fitMod
  browser()
}

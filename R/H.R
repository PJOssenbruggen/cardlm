#' \code{H} is a wrapper function for plotting a ideal bottleneck merge.
#'
#' @usage H(t0,t1,u1,u2,x1,x2) uses \code{xabparam}, \code{xab} and \code{uab} from the \code{cartools} Package.
# #' @examples
# #' H()
#' @export
H <- function(t0,t1,u1,u2,x1,x2) {
  dt <- t1 - t0
  Y  <- u2-u1
  Z  <- x2-x1-u1*dt
  A  <- dt
  B  <- dt^2/2
  C  <- -1*dt^2/2
  D  <- -dt^3/6
  a  <- (Y*D-Z*B)/(A*D-B*C)
  b  <- (Z*A-C*Y)/(A*D-B*C)
  df <- data.frame(t0,t1,u1,u2,x1,x2)
#  print(data.frame(a,b))
  browser()
  dlm()
  h <- dlm(FF = matrix(rep(1,length(df)), nrow = 1),
           V  = matrix(seq(1,length(df)), ncol = 1),
           GG = 1,
           W  = 2,
           m0 = rep(0,length(df)),
           C0 = 10 * diag(length(df))
           )
}

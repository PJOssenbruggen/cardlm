#' \code{cartools.example} is a wrapper function for estimating a successful, safe merge.
#'
#' @usage uses \code{xabparam}, \code{xab} and \code{uab} from the \code{cartools} Package.
# #' @examples
# #' cartools.example()
#' @export
cartools.example <- function() {
  lst     <- zipper.simulate(nveh1 = 5,nveh2 = 5,umn = 53.1,
            usd = 5, xstart1 = -700,delt = 0.125, tstart = 0,tend = 40,
            xfunnel = -500, leff = 14, size = 1,kfactor = 4/3)
  browser()
  df      <- as.matrix(lst[[6]])
  seq1    <- seq(1,dim(df)[1],10)
  u1      <- as.numeric(df[seq1,5])
  seq2    <- seq(2,dim(df)[1],10)
  u2      <- as.numeric(df[seq2,4])
  u1      <- ts(u1)

  deaths <- window(cartools[,2], end = 2007)
  population <- window(cartools[,6],end = 2007)
  Zt <- matrix(c(1,0),1,2)
  Ht <- matrix(NA)
  Tt <- matrix(c(1,0,1,1),2,2)
  Rt <- matrix(c(1,0),2,1)
  Qt <- matrix(NA)
  a1 <- matrix(c(1,0),2,1)
  P1 <- matrix(0,2,2)
  P1inf <- diag(2)
  model_gaussian <- SSModel(deaths/population ~ -1 + SSMcustom(Z = Zt, T=Tt,R=Rt,Q=Qt,a1=a1,P1=P1,P1inf=P1inf),
                            H = Ht)
  fit_guassian   <- fitSSM(model_gaussian, inits = c(0,0), method = "BFGS")
  out_guassian   <- KFS(fit_guassian$model)
  plot(deaths/population, col = gray(0.5))
  lines(out_guassian$a[,1], col = "red", lwd = 2)
  lines(out_guassian$muhat, col = "blue", lwd = 2)

}

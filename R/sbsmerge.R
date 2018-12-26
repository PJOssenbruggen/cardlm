#' \code{sbsmerge} is a wrapper function for estimating a successful, safe merge.
#'
#' @usage uses \code{xabparam}, \code{xab} and \code{uab} from the \code{cartools} Package.
# #' @examples
# #' sbsmerge()
#' @export
sbsmerge <- function() {
  nveh1   <- nveh2 <- 5
  umn     <- 53.1
  usd     <- 5
  xstart1 <- -700
  delt    <- 0.125
  tstart  <- 0
  tend    <- 40
  xfunnel <- -500
  leff    <- 14
  size    <- 20
  kfactor <- 4/3
  lst     <- zipper.simulate(nveh1,nveh2,umn,usd,xstart1,delt,tstart,tend,xfunnel,leff,size,kfactor)
  ### Establish data frames for vehicles 1 and 2, veh = 1 and veh = 2.
  df      <- as.matrix(lst[[6]])
  seq1    <- seq(1,dim(df)[1],10)
  u1      <- as.numeric(df[seq1,5])
  seq2    <- seq(2,dim(df)[1],10)
  u2      <- as.numeric(df[seq2,4])
  u1      <- ts(u1)
  u2      <- ts(u2)
  ### Establish a random walk filter model for veh = 1, mod1.
  build   <- function(parm) {dlm::dlmModPoly(1, dV = exp(parm[1]), dW = exp(parm[2]))}
  fit1    <- dlm::dlmMLE(u1, rep(0,2), build)
  parms1  <- unlist(build(fit1$par))[c("V","W")]
  dV1     <- parms1[1]
  dW1     <- parms1[2]
  rw1     <- dlm::dlmModPoly(order = 1, dV1, dW1)
  m0(rw1) <- umn*5280/3600
  C0(rw1) <- 60
  mod1    <- dlm::dlmFilter(u1, rw1)
  u1      <- c(umn*5280/3600, unlist(u1))
  df1     <- data.frame(u1 = ts(u1), m1 = mod1$m)
  se      <- rep(NA,length(size+1))
  k       <- seq(0,size)
  df1     <- data.frame(k, u1, m1 = as.numeric(mod1$m), se, m1up = se, m1lw = se)
  attach(mod1)
  for(i in 1:size+1) {
    df1$se[i] <- sqrt(unlist(dlm::dlmSvd2var(U.C[[i]], D.C[[i]])))
  }
  detach(mod1)
  df1$m1up <- df1$m1 + 1.96*df1$se
  df1$m1lw <- df1$m1 - 1.96*df1$se
  par(mfrow = c(2,2), pty = "s")
  plot(df1$k, df1$m1, xlab = "k", ylab = expression(u[k]), lwd = 2, col = "black", typ = "b",ylim = c(0,100))
  points(df1$k, df1$u1, col = gray(0.5), pch = 16)
  for(i in 1:length(k)) lines(c(df1$k[i], df1$k[i]), c(df1$m1lw[i], df1$m1up[i]))
  title(main = "Vehicle = 1")
  legend("topright", legend = c("Observed","Filtered", "95% C.I."),
  pch = c(16,1,1),
  lty = c(NA,NA,1),
  col = c(gray(0.5),"black",gray(0.0)),
  bty = "n"
  )
  ### Establish a random walk filter model for veh = 2.
  fit2    <- dlm::dlmMLE(u2, rep(0,2), build)
  parms2  <- unlist(build(fit2$par))[c("V","W")]
  dV2     <- parms2[1]
  dW2     <- parms2[2]
  rw2     <- dlm::dlmModPoly(order = 1, dV2, dW2)
  m0(rw2) <- umn*5280/3600
  C0(rw2) <- 60
  mod2    <- dlm::dlmFilter(u2, rw2)
  u2      <- c(umn*5280/3600, unlist(u2))
  df2     <- data.frame(u2 = ts(u2), m2 = mod2$m)
  se      <- rep(NA,length(size+1))
  k       <- seq(0,size)
  df2     <- data.frame(k, u2, m2 = as.numeric(mod2$m), se, m2up = se, m2lw = se)
  attach(mod2)
  for(i in 1:size+1) {
    df2$se[i] <- sqrt(unlist(dlm::dlmSvd2var(U.C[[i]], D.C[[i]])))
  }
  detach(mod2)
  df2$m2up <- df2$m2 + 1.96*df2$se
  df2$m2lw <- df2$m2 - 1.96*df2$se
  plot(df2$k, df2$m2, xlab = "k", ylab = expression(u[k]), lwd = 2, col = "black", typ = "b",ylim = c(0,100))
  points(df2$k, df2$u2, col = gray(0.5), pch = 16)
  for(i in 1:length(k)) lines(c(df2$k[i], df2$k[i]), c(df2$m2lw[i], df2$m2up[i]))
  title(main = "Vehicle = 2")
  ### Plot results
  df12    <- cbind(df1,df2)
  u1      <- seq(40,80,0.125)
  mn1     <- df12[21,3]
  sd1     <- df12[21,4]
  mn2     <- df12[21,9]
  sd2     <- df12[21,10]
  u1plot  <- dnorm(u1,mn1,sd1)
  u2plot  <- dnorm(u1,mn2,sd2)
  browser()
  plot(u1,u2plot, lwd = 2, typ = "l", col = "red", xlab = "u, speed (fps)", ylab = "", ylim = c(0,0.2))
  lines(u1,u1plot,col = "black", lwd = 2)
  return(df12)
}

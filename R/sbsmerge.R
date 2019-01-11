#' \code{sbsmerge} is a wrapper function for estimating a successful, safe merge.
#'
#' @usage uses \code{xabparam}, \code{xab} and \code{uab} from the \code{cartools} Package.
# #' @examples
# #' sbsmerge()
#' @export
sbsmerge <- function() {
  set.seed(123)
  nveh1   <- nveh2 <- 5
  umn     <- 53.1
  usd     <- 5
  xstart1 <- -700
  delt    <- 0.125
  tstart  <- 0
  tend    <- 40
  xfunnel <- -500
  leff    <- 14
  size    <- 100
  kfactor <- 4/3
  lst     <- sbs.simulate(nveh1,nveh2,umn,usd,xstart1,delt,tstart,tend,xfunnel,leff,size,kfactor)
  ### Establish data frames for vehicles 1 and 2, veh = 1 and veh = 2 where u1 and u2 are the speeds of veh 1 and 2.
  ### Dt = t2 - t1
### Lead Vehicle 1, Following vehicle 2 ##########################################
  df      <- as.matrix(lst[[6]])
  head(df)
  for(i in 1:10) {
    stp   <- seq(i,dim(df)[1],10)
    te    <- as.numeric(df[stp,2])
    ue    <- as.numeric(df[stp,5])
    td    <- as.numeric(df[stp,3])
    ud    <- as.numeric(df[stp,4])
    tmp   <- data.frame(te,ue,td,ud)
    if(i == 1) df.reg <- tmp else df.reg <- cbind(df.reg, tmp)
  }
  hdr <- rep(seq(1,10),4)
  o   <- order(hdr)
  hdr <- hdr[o]
  names(df.reg) <- paste0(names(df.reg),hdr)
  regress <- function(lead.veh, df.reg) {
    cols  <- matrix(seq(1,40),10,4,byrow = TRUE)
    par(mfrow = c(1,2), pty = "s")
    data <- df.reg[,cols[lead.veh + 1,]]
    data <- data.frame(dt = data[,3] - data[,1], ue = data[,2], u0 = data[,4])
    plot(data$ue, data$dt,
         ylab = expression(Delta*"t"), xlab = expression(u[e]),
         ylim = c(0,1.1*max(data$dt)),
         xlim = c(min(data$ue), max(data$ue))
    )
    title(main = "Predictions", sub = paste("Leading Vehicle:",lead.veh))
    dt.u.reg  <- lm(dt ~ ue, data)
    summary(dt.u.reg)
    lm.prd  <- predict(dt.u.reg,  data, interval = "prediction")
    lines(data$ue, as.numeric(lm.prd[,1]), lwd = 2, col = "red")
    lines(data$ue, as.numeric(lm.prd[,2]), col = "red", lty = 2)
    lines(data$ue, as.numeric(lm.prd[,3]), col = "red", lty = 2)
    plot(data$ue, data$u0,
         ylab = expression(u[0]), xlab = expression(u[e]),
         ylim = c(0,1.1*max(data$u0)),
         xlim = c(min(data$ue), max(data$ue))
    )
    u0.u.reg  <- lm(u0 ~ ue, data)
    summary(u0.u.reg)
    lm.prd  <- predict(u0.u.reg,  data, interval = "prediction")
    lines(data$ue, as.numeric(lm.prd[,1]), lwd = 2, col = "blue")
    lines(data$ue, as.numeric(lm.prd[,2]), col = "blue", lty = 2)
    lines(data$ue, as.numeric(lm.prd[,3]), col = "blue", lty = 2)
    title(main = "Predictions", sub = paste("Following Vehicle:",lead.veh + 1))
  }
  regress(lead.veh = 2, df.reg)

  data("alcohol")
  deaths <- window(alcohol[,2], end = 2007)
  population <- window(alcohol[,6], end = 2007)
  dt <- 0.125
  Zt <- matrix(c(1,0), 1, 2)
  Ht <- matrix(NA)
  Tt <- matrix(c(1,0,1,1),2,2)
  Rt <- matrix(c(1,0),2,1)
  Qt <- matrix(NA)
  a1 <- matrix(c(1,0),2,1)
  P1 <- matrix(0,2,2)
  P1inf <- diag(2)
  par(mfrow = c(1,2), pty = "s")
  plot(deaths/population, col = gray(0.5), typ = "b", ylim = c(0,80))
  model_gaussian <- SSModel(deaths/population ~ -1 + SSMcustom(Z = Zt, T = Tt,
              R = Rt, Q = Qt, a1 = a1, P1 = P1, P1inf = P1inf), H = Ht)
  fit_gaussian   <- fitSSM(model_gaussian, inits = c(0,0), method = "BFGS")
  out_gaussian   <- KFS(fit_gaussian$model)
  print(out_gaussian)
  attributes(out_gaussian)
  lines(out_gaussian$att[,1], col = "red", lwd = 2)
  lines(out_gaussian$muhat, col = "blue", lwd = 2)
  title(main = "Alcohol", sub = "Gaussian")
  legend("topleft", legend = c("Observations", "One-step head pred.", "Smoothed estimates"),
         lty = c(NA,1,1), col = c(gray(0.5),"red","blue"), lwd = c(NA,2,2), pch = c(1,NA,NA), bty = "n"
  )
  P <- ts(as.matrix(out_gaussian$P), frequency = 1, start = 1969, end = 2007)
  plot(P, col = "black", typ = "l", lwd = 2)
  title("Error Variance")
  browser()


  ### dlmModReg for u2 ~ dt
  if(FALSE) {
    build.reg<- function(parm) {dlm::dlmModReg(Dt.ts, dV = exp(parm[1]), dW = exp(parm[2:3]))}
    outMLE   <- dlmMLE(u2, parm = rep(0,3), build.reg)
    exp(outMLE$par)
    outMLE$value
    mod.reg <- build.reg(outMLE$par)
    outS    <- dlm::dlmSmooth(u2, mod.reg)
    plot(dropFirst(outS$s))
    outF    <- dlm::dlmFilter(u2, mod.reg)
    plot(outF$m)
  }
  ### dlmModReg for dt ~ u1
  if(FALSE) {
    build.reg<- function(parm) {dlm::dlmModReg(u1, dV = exp(parm[1]), dW = exp(parm[2:3]))}
    outMLE   <- dlmMLE(Dt.ts, parm = rep(0,3), build.reg)
    exp(outMLE$par)
    outMLE$value
    mod.reg <- build.reg(outMLE$par)
    outS    <- dlm::dlmSmooth(Dt.ts, mod.reg)
    plot(dropFirst(outS$s))
    outF    <- dlm::dlmFilter(Dt.ts, mod.reg)
    plot(dropFirst(outS$s))
  }
  ##########################################################################################
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
#  par(mfrow = c(2,2), pty = "s")
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
  browser()
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

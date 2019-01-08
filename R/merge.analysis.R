#' \code{merge.analysis} is a wrapper function for estimating a successful, safe merge.
#'
#' @usage uses \code{merge.analysis}, \code{xab} and \code{uab} from the \code{cartools} Package.
#' @param tdf1df2, time, speed and locations, a matrix
# #' @examples
# #' merge.analysis(tdf1df2)
#' @export
merge.analysis <- function(tdf1df2) {
  par(mfrow = c(1,1), pty = "s")
  plot(tdf1df2[,1], tdf1df2[,3], type = "l", ylim = c(-1200,2400),
       ylab = expression("x, feet"), xlab = "t, seconds")
  xveh  <- seq(3,30,3)
  uveh  <- seq(2,29,3)
  tveh  <- tdf1df2[,1]
  h     <- x.xe <- u.x0  <- u.xe <- t.x0 <- t.xe <- rep(NA,10)
  u     <- tdf1df2[tdf1df2[,xveh[1]] <= 0, uveh[1]]
  t     <- tveh[tdf1df2[,xveh[1]] <= 0]
  u.x0[1] <- u[length(t)]
  t.x0[1] <- t[length(t)]
  x     <- tdf1df2[tdf1df2[,xveh[1]] <= -500, xveh[1]]
  u     <- tdf1df2[tdf1df2[,xveh[1]] <= -500, uveh[1]]
  t     <- tveh[tdf1df2[,xveh[1]] <= -500]
  x.xe[1] <- x[length(t)]
  u.xe[1] <- u[length(t)]
  t.xe[1] <- t[length(t)]
  for(i in 2:10) {
    lines(tdf1df2[,1], tdf1df2[,xveh[i]])
    u     <- tdf1df2[tdf1df2[,xveh[i]] <= 0, uveh[i]]
    t     <- tveh[tdf1df2[,xveh[i]] <= 0]
    u.x0[i] <- u[length(t)]
    t.x0[i] <- t[length(t)]
    x     <- tdf1df2[tdf1df2[,xveh[i]] <= -500, xveh[i]]
    u     <- tdf1df2[tdf1df2[,xveh[i]] <= -500, uveh[i]]
    t     <- tveh[tdf1df2[,xveh[i]] <= -500]
    x.xe[i] <- x[length(t)]
    u.xe[i] <- u[length(t)]
    t.xe[i] <- t[length(t)]
  }
  for(i in 1:10) h[i] <- hsafe(u.x0[i], 14)
  abline(h = c(0,-500), col = gray(0.5))
  abline(v = 0, col = gray(0.5))
  axis(side = 4, at = 0, label = expression(x[0]))
  axis(side = 4, at = -500, label = expression(x[e]))
  df <- data.frame(t.xe, u.xe, t.x0, u.x0, h, x.xe)
  for(i in 1:10) points(t.x0[i], -h[i])
  dt.0 <- t.x0 - t.xe
  df <- cbind(df, dt.0)
  for(i in 2:10) {
    lines(c(t.xe[i],t.x0[i]), c(x.xe[i], x.xe[i] + u.xe[i] * df[i,7]), lwd = 3, col = "red")
  }
  title(main = "SBS Driver Inefficiencies", sub = expression("Decision location: "*x[e]))

### u.1 #####################################################################
  browser()
  u    <- ts(data.frame(u1 =  tdf1df2[,2], u2 =  tdf1df2[,5],
                           x1 =  tdf1df2[,3], x2 =  tdf1df2[,6]),
                start = 0, end = 40, frequency = 8)
  dt <- 0.125
  Zt <- matrix(1)
  Ht <- matrix(NA)
  Tt <- matrix(1)
  Rt <- matrix(1)
  Qt <- matrix(NA)
  a1 <- matrix(1)
  P1 <- matrix(0)
  P1inf <- diag(1)
  start <- c(0,1)
  end   <- c(36,0)
  uPred <- window(u, start, end)
  par(mfrow = c(1,1), pty = "s")
  model_u <- SSModel(u1 ~ -1 +
                SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, a1 = a1, P1 = P1, P1inf = P1inf),
                H = Ht, data = uPred)
  fit_u   <- fitSSM(model_u, inits = c(0,0), method = "BFGS")
  out_u   <- KFS(fit_u$model)
  print(out_u)
  pred <- predict(fit_u$model,
            newdata = SSModel(ts(matrix(NA,nprd,1), start = end) ~ -1 +
            SSMcustom(Z = fit_u$model$Z, T = fit_u$model$T, R = fit_u$model$R, Q = fit_u$model$Q),
            distribution = "gaussian"),
              interval = "confidence", nsim = 10000)
  plot(u[,1],  col = gray(0.5), type = "p",
       ylab = "u, feet per second", xlab = "t, seconds",
       lwd = 2, ylim = c(min(u[,1] - 20), max(u[,1])))
  lines(out_u$att, col = "red", lwd = 4)
  lines(out_u$muhat, col = "yellow", lwd = 2)
  lines(pred[,1], col = "black", lwd = 3, lty = 1)
  lines(pred[,2], col = "black", lwd = 2, lty = 3)
  lines(pred[,3], col = "black", lwd = 2, lty = 3)
  title(main = "Lead Vehicle", sub = "Gaussian")
  legend("bottomleft",
         legend = c(
           expression(u[1]),
           "Filtered estimates", "Smoothed estimates", "95% C.I."
         ),
         pch = c(1,NA,NA,NA,NA,NA),
         lty = c(NA,1,1,1,3),
         lwd = c(NA,4,2,3,2),
         col = c("blue","red", "yellow", "black", "black"),
         bty = "n")

  browser()

### model_u1rw #####################################################################
  dt    <- 0.125
  Zt    <- matrix(c(1,0),1,2)
  Ht    <- matrix(NA)
  Tt    <- matrix(c(1,0,1,1),2,2)
  Rt    <- matrix(c(1,0),2,1)
  Qt    <- matrix(NA)
  a1    <- matrix(c(1,0),2,1)
  P1    <- matrix(0,2,2)
  P1inf <- diag(2)
  start <- c(0,1)
  end   <- c(30,0)
  nprd  <- 80
  uPred <- window(u, start, end)
  model_u1rw <- SSModel(u1 ~ -1 +
             SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, a1 = a1, P1 = P1, P1inf = P1inf),
             H = Ht, data = uPred)
  fit_urw   <- fitSSM(model_u1rw, inits = c(0,0), method = "BFGS")
  out_urw   <- KFS(fit_urw$model)
  predrw    <- predict(fit_urw$model,
          newdata = SSModel(ts(matrix(NA,nprd,1), start = end) ~ -1 +
          SSMcustom(Z = fit_urw$model$Z, T = fit_urw$model$T,
                    R = fit_urw$model$R, Q = fit_urw$model$Q),
          distribution = "gaussian"),
          interval = "confidence", nsim = 10000)
  par(mfrow = c(1,1), pty = "s")
  plot(u[,1],  col = gray(0.5), type = "p",
       ylab = "u, feet per second", xlab = "t, seconds",
       lwd = 2, ylim = c(min(u[,1] - 20), max(u[,1])))
  lines(out_urw$att[,1], col = "red", lwd = 4)
  lines(out_urw$muhat, col = "yellow", lwd = 2)
  lines(predrw[,1], col = "black", lwd = 3, lty = 1)
  lines(predrw[,2], col = "black", lwd = 2, lty = 3)
  lines(predrw[,3], col = "black", lwd = 2, lty = 3)
  title(main = "Lead Vehicle Random Walk", sub = "Gaussian")
  legend("bottomleft",
         legend = c(
           expression(u[1]),
           "Filtered estimates", "Smoothed estimates", "95% C.I."
         ),
         pch = c(1,NA,NA,NA,NA,NA),
         lty = c(NA,1,1,1,3),
         lwd = c(NA,4,2,3,2),
         col = c(gray(0.5),"red", "yellow", "black", "black"),
         bty = "n")
  browser()
  out_urw$a
  out_urw$Pinf
  out_urw$att

### model_x1rw #####################################################################
  dt    <- 0.125
  Zt    <- matrix(c(0,1,0),1,3)
  Ht    <- matrix(NA)
  Tt    <- matrix(c(1,0,0,dt,1,0,0,1,1),3,3)
  Rt    <- matrix(c(0,1,0),3,1)
  Qt    <- matrix(NA)
  a1    <- matrix(c(0,1,0),3,1)
  P1    <- matrix(0,3,3)
  P1inf <- diag(3)
  start <- c(0,1)
  end   <- c(30,0)
  uPred <- window(u, start, end)
  par(mfrow = c(1,1), pty = "s")
  model_x1rw <- SSModel(x1 ~ -1 +
                  SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, a1 = a1, P1 = P1, P1inf = P1inf),
                      H = Ht, data = uPred)
  ts.plot()
  fit_x1rw   <- fitSSM(model_x1rw, inits = c(0,0,0), method = "BFGS")
  out_x1rw   <- KFS(fit_x1rw$model)
  print(out_x1rw)
  print(attributes(out_x1rw))
  predrw <- predict(fit_x1rw$model,
                    newdata = SSModel(ts(matrix(NA,nprd,1), start = end) ~ -1 +
                                SSMcustom(Z = fit_x1rw$model$Z, T = fit_x1rw$model$T,
                                R = fit_x1rw$model$R, Q = fit_x1rw$model$Q),
                                distribution = "gaussian"),
                    interval = "confidence", nsim = 10000)
  plot(u[,3],  col = gray(0.5), type = "p",
       ylab = "u, feet per second", xlab = "t, seconds",
       lwd = 2)
  lines(out_x1rw$att[,2], col = "red", lwd = 4)
  lines(out_x1rw$muhat, col = "yellow", lwd = 2)
  lines(predrw[,1], col = "black", lwd = 3, lty = 1)
  lines(predrw[,2], col = "black", lwd = 2, lty = 3)
  lines(predrw[,3], col = "black", lwd = 2, lty = 3)
  title(main = "Lead Vehicle Random Walk 2", sub = "Gaussian")
  legend("bottomleft",
         legend = c(
           expression(u[1]),
           "Filtered estimates", "Smoothed estimates", "95% C.I."
         ),
         pch = c(1,NA,NA,NA,NA,NA),
         lty = c(NA,1,1,1,3),
         lwd = c(NA,4,2,3,2),
         col = c("blue","red", "yellow", "black", "black"),
         bty = "n")
  browser()

  ### model_x2rw #####################################################################
  dt    <- 0.125
  Zt    <- matrix(c(0,1,0),1,3)
  Ht    <- matrix(NA)
  Tt    <- matrix(c(1,0,0,dt,1,0,0,1,1),3,3)
  Rt    <- matrix(c(0,1,0),3,1)
  Qt    <- matrix(NA)
  a1    <- matrix(c(0,1,0),3,1)
  P1    <- matrix(0,3,3)
  P1inf <- diag(3)
  start <- c(0,1)
  end   <- c(36,0)
  uPred <- window(u, start, end)
  par(mfrow = c(1,1), pty = "s")
  model_x2rw <- SSModel(u1 ~ -1 +
                SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, a1 = a1, P1 = P1, P1inf = P1inf),
                H = Ht, data = uPred)
  fit_x2rw   <- fitSSM(model_x2rw, inits = c(0,0,0), method = "BFGS")
  out_x2rw   <- KFS(fit_x2rw$model)
  print(out_x2rw)
  print(attributes(out_x2rw))
  predrw <- predict(fit_x2rw$model,
                    newdata = SSModel(ts(matrix(NA,nprd,1), start = end) ~ -1 +
                                SSMcustom(Z = fit_x2rw$model$Z, T = fit_x2rw$model$T,
                                R = fit_x2rw$model$R, Q = fit_x2rw$model$Q),
                                distribution = "gaussian"),
                    interval = "confidence", nsim = 10000)
  tseq  <- seq(tsp(u)[1], tsp(u)[2], dt)
  tseq1 <- seq(start, end, dt)
  tseq2 <- seq(end+1, end + nprd, dt)
  par(mfrow = c(1,2), pty = "s")
  plot(tseq, as.numeric(u[,1]), type = "n", col = gray(0.5),
       ylab = expression(u[1]*"(t), feet per second"), xlab = "t, seconds",
       lwd = 2, ylim = c(min(u[,1] - 20), max(u[,1])))
  points(tseq, as.numeric(u[,1]), col = "blue")
  lines(out_x2rw$att[,2], col = "red", lwd = 4)
  lines(out_x2rw$muhat, col = "yellow", lwd = 2)
  lines(tseq2, as.numeric(predrw[,1]), col = "black", lwd = 3, lty = 1)
  lines(tseq2, as.numeric(predrw[,2]), col = "black", lwd = 2, lty = 3)
  lines(tseq2, as.numeric(predrw[,3]), col = "black", lwd = 2, lty = 3)
  title(main = "Lead Vehicle Random Walk 2", sub = "Gaussian")
  legend("bottomleft",
         legend = c(
           expression(u[1]),
           "Filtered estimates", "Smoothed estimates", "95% C.I."
         ),
         pch = c(1,NA,NA,NA,NA,NA),
         lty = c(NA,1,1,1,3),
         lwd = c(NA,4,2,3,2),
         col = c("blue","red", "yellow", "black", "black"),
         bty = "n")

  plot(tseq, as.numeric(u[,3]), type = "l", col = gray(0.5),
       ylab = expression(hat(x)[1]*"(t), feet"), xlab = "t, seconds",
       lwd = 2, ylim = c(-800, max(u[,3])))
  lines(out_x2rw$att[,1]-rep(700,length(tseq)), col = "red", lwd = 4)
  abline(h = c(0,-500), col = gray(0.5))
  abline(v = 0, col = gray(0.5))
  title(main = "No Random Effects for Location", sub = "Gaussian")

  browser()
### model_temp ######################################################################
  # Example of multivariate local level model with only one state
  # Two series of average global temperature deviations for years 1880-1987
  # See Shumway and Stoffer (2006), p. 327 for details
  par(mfrow = c(1,1), pty = "s")
  data("GlobalTemp")
  model_temp <- SSModel(GlobalTemp ~ SSMtrend(1, Q = NA, type = "common"),
                        H = matrix(NA, 2, 2))
  # Estimating the variance parameters
  inits      <- chol(cov(GlobalTemp))[c(1, 4, 3)]
  inits[1:2] <- log(inits[1:2])
  fit_temp   <- fitSSM(model_temp, c(0.5*log(.1), inits), method = "BFGS")
  out_temp   <- KFS(fit_temp$model)
  ts.plot(cbind(model_temp$y, coef(out_temp)), col = c("red","blue","black"),
          lwd = c(1,1,4))
  legend("bottomright",
         legend = c(colnames(GlobalTemp), "Smoothed signal"),
         col = c("red","blue","black"), lty = 1, lwd = c(1,1,4))
  title(main = "global temperature", sub = "observations from two series.")
  browser()

### model_u12 ##########################################################
  u12        <- ts(data.frame(u1 =  tdf1df2[,2], u2 =  tdf1df2[,5]),
                   start = 0, end = 40, frequency = 8)
  model_u12  <- SSModel(u12 ~ SSMtrend(1, Q = NA, type = "common"),
                        H = matrix(NA, 2, 2), data = u12)
  inits      <- chol(cov(u12))[c(1, 4, 3)]
  inits[1:2] <- log(inits[1:2])
  fit_u12    <- fitSSM(model_u12, c(0.5*log(.1), inits), method = "BFGS")
  out_u12    <- KFS(fit_u12$model)
  ts.plot(cbind(model_u12$y, coef(out_u12)), col = c("red","blue","black"),
          lwd = c(1,1,4), ylab = expression(hat(u)[2](t)))
  legend("topright",
         legend = c(expression(u[1]), expression(u[2]), "Smoothed signal"),
         col = c("red","blue","black"), lty = 1, lwd = c(1,1,4))
  title(main = "Trend", sub = "observations from two series.")
  browser()

### model_u12v1 ###########################################################
  u12<- as.ts(data.frame(u1 =  tdf1df2[,2], u2 =  tdf1df2[,5]),
              start = 0, end = 40, frequency = 8)
  Zt <- matrix(c(1,1),2,2)
  Ht <- matrix(c(NA,NA), 2, 2)
  Tt <- matrix(c(1,0,1,1), 2, 2)
  Rt <- matrix(c(1, 0), 2, 1)
  Qt <- matrix(NA)
  a1 <- matrix(c(0, 0), 2, 1)
  P1 <- matrix(0,2,2)
  P1inf <- diag(2)
  start <- 1
  nprd  <- 25
  end   <- 321 - nprd
  uPred <- window(u12, start, end)
  model_u12v1 <- SSModel(uPred ~ -1 +
                SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, a1 = a1, P1 = P1, P1inf = P1inf),
                H = Ht, data = uPred)
  print(model_u12v1)
  inits      <- chol(cov(u12))[c(1, 4, 3)]
  inits[1:2] <- log(inits[1:2])
  fit_u12v1  <- fitSSM(model_u12v1, c(0.5*log(.1), inits), method = "BFGS")
  out_u12v1  <- KFS(fit_u12v1$model)
  ts.plot(cbind(model_u12v1$y, coef(out_u12v1))[,1:3], col = c("red","blue","black"),
          lwd = c(1,1,4), ylim = c(40,120), ylab = expression(hat(u)[2](t)))
  legend("topright",
         legend = c(expression(u[1]), expression(u[2]), "Smoothed signal"),
         col = c("red","blue","black"), lty = 1, lwd = c(1,1,4))
  title(main = "Durbin/Koopman model", sub = "observations from two series.")
  browser()

### model_u3 #####################################################################
  u    <- as.ts(data.frame(u1 =  tdf1df2[,2], u2 =  tdf1df2[,5],
                           x1 =  tdf1df2[,3], x2 =  tdf1df2[,6]))
  dt <- 0.125
  Zt <- matrix(c(1,0), 1, 2)
  Ht <- matrix(NA)
  Tt <- matrix(c(1,0,1,1),2,2)
  Rt <- matrix(c(1,0),2,1)
  Qt <- matrix(NA)
  a1 <- matrix(c(1,0),2,1)
  P1 <- matrix(0,2,2)
  P1inf <- diag(2)
  par(mfrow = c(1,1), pty = "s")
  model_u3 <- SSModel(u1 ~ -1 +
                SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, a1 = a1, P1 = P1, P1inf = P1inf),
                H = Ht, data = u)
  fit_u   <- fitSSM(model_u3, inits = c(0,0), method = "BFGS")
  out_u   <- KFS(fit_u$model)
  print(out_u)
  print(attributes(out_u))

  plot(u[,1], ylim = c(min(c(u[,1],u[,2])) -20, max(c(u[,1],u[,2]))),
       col = "blue", type = "p", ylab = expression(hat(u)[1](t)))
  points(u[,2], col = "orange")
  lines(out_u$att[,1], col = "red", lwd = 4)
  lines(out_u$muhat, col = "yellow", lwd = 2)
  title(main = "Lead Vehicle", sub = "Gaussian: one series")
  legend("bottomleft", legend = c("Lead", "Following", "Filtered estimates","Smoothed estimates"),
         lty = c(NA,NA,1,1),
         col = c("blue", "orange", "red", "yellow"),
         lwd = c(NA,NA,4,2),
         pch = c(1,1,NA,NA), bty = "n")
  browser()

### model_u2 ########################################################################
  start <- 1
  nprd  <- 25
  end   <- 321 - nprd
  uPred <- window(u, start, end)
  par(mfrow = c(1,1), pty = "s")
  ts.plot(uPred[,3:4], lty = c(1,1), lwd = c(2,2), col = c("blue","red"), ylab = "x, feet")
  legend("topleft",  legend = colnames(uPred)[3:4],
         lty= c(1,1), lwd = c(2,2), col = c("blue","red"),
         bty = "n")
  browser()
  model_u2  <- SSModel(uPred[,3:4] ~
                  SSMtrend(2, Q = list(matrix(NA,2,2), matrix(0,2,2))) +
                  SSMcustom(Z = diag(1,2), T = diag(0,2), Q = matrix(NA,2,2),
                                           P1 = matrix(NA,2,2)),
                  distribution = "gaussian"
                  )
  updatefn <- function(pars, model, ...) {
    Q <- diag(pars[1:2])
    Q[upper.tri(Q)] <- pars[3]
    model["Q", etas = "level"] <- crossprod(Q)
    Q <- diag(pars[4:5])
    Q[upper.tri(Q)] <- pars[6]
    model["Q", etas = "custom"] <- model["P1", states = "custom"] <- crossprod(Q)
    model
  }
  init     <- chol(cov(u[,3:4]))
  fitinit  <- fitSSM(model_u2, updatefn = updatefn,
                     inits = rep(c(diag(init), init[upper.tri(init)]),2),
                     method = "BFGS")
  print(-fitinit$optim.out$val)
  fit <- fitSSM(model_u2, updatefn = updatefn,
                inits = fitinit$optim.out$par,
                method = "BFGS", nsim = 250)
  print(-fitinit$optim.out$val)
  varcor <- fit$model["Q", etas = "level"]
  varcor[upper.tri(varcor)] <- cov2cor(varcor)[upper.tri(varcor)]
  print(varcor, digits = 2)
  varcor <- fit$model["Q", etas = "custom"]
  varcor[upper.tri(varcor)] <- cov2cor(varcor)[upper.tri(varcor)]
  print(varcor, digits = 2)
  out <- KFS(fit$model, nsim = 1000)
  print(out)
  plot(coef(out, states = c("level", "custom")), main = "Smoothed states", yax.flip = TRUE)
  browser()
  res <- rstandard(KFS(fit$model))
  acf(res, na.action = na.pass)
  pred <- predict(fit$model,
            newdata = SSModel(ts(matrix(NA,nprd,2), start = end) ~ -1 +
                SSMcustom(Z = fit$model$Z, T = fit$model$T, R = fit$model$R, Q = fit$model$Q),
                distribution = "gaussian"),
                interval = "confidence", nsim = 10000)
  trend <- signal(out, "trend")$signal
  browser()
  tseq  <- seq(tsp(u)[1], tsp(u)[2], dt)
  tseq1 <- seq(start, end, dt)
  tseq2 <- seq(end+1, end + nprd, dt)

  plot(tseq, as.numeric(u[,3]), type = "n", col = gray(0.5),
       ylab = "x, feet", xlab = "t, seconds",
       lwd = 2, ylim = c(-700, max(c(u[,3],u[,4]))),
       main = "Car Following")
  for(i in 1:2) {
    if(i == 1) {
      points(tseq, as.numeric(u[,3]), col = "blue")
      lines(tseq1, as.numeric(trend[,1]), lwd = 2, col = "blue")
      lines(tseq2, as.numeric(pred[[1]][,1]), col = "black", lwd = 3, lty = 1)
      lines(tseq2, as.numeric(pred[[1]][,2]), col = "black", lwd = 2, lty = 3)
      lines(tseq2, as.numeric(pred[[1]][,3]), col = "black", lwd = 2, lty = 3)
    } else {
      points(tseq, as.numeric(u[,4]), col = "red")
      lines(tseq1, as.numeric(trend[,2]), lwd = 2, col = "red")
      lines(tseq2, as.numeric(pred[[2]][,1]), col = "black", lwd = 3, lty = 1)
      lines(tseq2, as.numeric(pred[[2]][,2]), col = "black", lwd = 2, lty = 3)
      lines(tseq2, as.numeric(pred[[2]][,3]), col = "black", lwd = 2, lty = 3)
    }
    abline(v = 0, col = gray(0.5))
    abline(h = 0, col = gray(0.5))
    legend("topleft",
           legend = c(
             expression(x[1]),
             expression(x[2]),
             "filtered", "predictions", "95% C.I."
           ),
           pch = c(1,1,NA,NA,NA),
           lty = c(NA,NA,1,1,3),
           lwd = c(NA,NA,2,3,2),
           col = c("blue","red", "black", "black", "black"),
           bty = "n")
  }

#########################################################################
}


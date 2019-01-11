#' \code{cartools.example} is a wrapper function for estimating a successful, safe merge.
#'
#' @usage uses \code{xabparam}, \code{xab} and \code{uab} from the \code{cartools} Package.
# #' @examples
# #' cartools.example(u,x)
#' @export
cartools.example <- function(u,x) {

  ### Gyro y_t = Z_t*alpha_t + epsilon; epsilon ~ N(0, sigma_epsilon^2)
  ### alpha = T_t*mu + R_t*nu, Q ~ N(0, sigma_nu^2)
  set.seed(123)
  obs <- BM(x = 1, t0 = 0, T = 3, N = 12)
  alpha <- -1
  sigma_epsilon   <- 5
  Q   <- sigma_nu <- 6
  Zt  <- matrix(alpha)
  Ht  <- matrix(sigma_epsilon)
  Tt  <- matrix(1)
  Rt  <- matrix(1)
  Qt  <- matrix(sigma_nu)
  a1  <- matrix(0)
  P1  <- matrix(0)
  P1inf <- matrix(1)
  par(mfrow = c(1,2), pty = "s")
  plot(obs, col = gray(0.5), typ = "b", ylim = c(-2,3))
  model_gyro <- SSModel(obs ~ -1 + SSMcustom(Z = Zt, T = Tt,
               R = Rt, Q = Qt, a1 = a1, P1 = P1, P1inf = P1inf), H = Ht)
  attributes(model_gyro)
  if(is.na(Q)) {
    fit_gyro   <- fitSSM(model_gyro, inits = c(0,0), method = "BFGS")
    out_gyro   <- KFS(fit_gyro$model)
  } else {
    out_gyro   <- KFS(model_gyro)
    print(model_gyro$Q)
    print(model_gyro$H)
    }
  attributes(out_gyro)
  print(out_gyro)
  lines(out_gyro$att, col = "red", lwd = 2)
  lines(out_gyro$muhat, col = "blue", lwd = 2)
  title(main = "Gyro", sub = "Gaussian")
  legend("topleft", legend = c("Observations", "One-step head pred.", "Smoothed estimates"),
         lty = c(NA,1,1), col = c(gray(0.5),"red","blue"), lwd = c(NA,2,2), pch = c(1,NA,NA), bty = "n"
  )
  P <- ts(as.matrix(out_gyro$P), frequency = 4, start = 0, end = 3)
  plot(P, col = "black", typ = "l", lwd = 2)
  title("Error Variance")
  browser()


###############################################################################

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

############################################################
  data("alcohol")
  start <- 1969
  end   <- 2007
  alcoholPred          <- window(alcohol, start, end)
  par(mfrow = c(1,1), pty = "s")
  ts.plot(alcohol[,1:2], lty = c(1,1), lwd = c(2,2), col = c("blue","red"), ylab = "deaths")
  legend("topleft",  legend = colnames(alcohol)[1:2],
         lty= c(1,1), lwd = c(2,2), col = c("blue","red"),
         bty = "n")
  model_alc_gauss2  <- SSModel(alcoholPred[,1:2] ~
                              SSMtrend(2, Q = list(matrix(NA,2,2), matrix(0,2,2))) +
                              SSMcustom(Z = diag(1,2), T = diag(0,2), Q = matrix(NA,2,2),
                                           P1 = matrix(NA,2,2)),
                              distribution = "gaussian",
                              u = alcoholPred[,5:6])
  browser()
  updatefn <- function(pars, model, ...) {
    Q <- diag(pars[1:2])
    Q[upper.tri(Q)] <- pars[3]
    model["Q", etas = "level"] <- crossprod(Q)
    Q <- diag(pars[4:5])
    Q[upper.tri(Q)] <- pars[6]
    model["Q", etas = "custom"] <- model["P1", states = "custom"] <- crossprod(Q)
    model
  }
  init     <- chol(cov(alcoholPred[,1:2]))
  fitinit  <- fitSSM(model_alc_gauss2, updatefn = updatefn,
                     inits = rep(c(diag(init), init[upper.tri(init)]),2),
                     method = "BFGS")
  print(-fitinit$optim.out$val)
  fit <- fitSSM(model_alc_gauss2, updatefn = updatefn,
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
  plot(coef(out, states = c("level", "custom")), main = "Smoothed states",yax.flip = TRUE)
  res <- rstandard(KFS(fit$model))
  acf(res, na.action = na.pass)
  pred <- predict(fit$model,
                  newdata <- SSModel(ts(matrix(NA,5,2), start = end) ~ -1 +
                            SSMcustom(Z = fit$model$Z, T = fit$model$T, R = fit$model$R, Q = fit$model$Q),
                                     u = 1, distribution = "gaussian"),
                  interval = "confidence", nsim = 10000)
  trend <- signal(out, "trend")$signal
  par(mfrow = c(1,1), pty = "s")
  tseq  <- seq(tsp(alcohol)[1], tsp(alcohol)[2])
  tseq1 <- seq(start,end)
  tseq2 <- seq(end+1, end + 5)
  plot(tseq, as.numeric(alcohol[,1]), typ = "n", col = gray(0.5),
       xlab = NULL, ylab = "deaths", lwd = 2, ylim = c(0, 700),
       main = "Alcohol Deaths in Finland")
  for(i in 1:2) {
    if(i == 1) {
      points(tseq, as.numeric(alcohol[,1]), col = "blue")
      lines(tseq1, as.numeric(trend[,1]), lwd = 2, col = "blue")
      lines(tseq2, as.numeric(pred[[1]][,1]), col = "black", lwd = 3, lty = 1)
      lines(tseq2, as.numeric(pred[[1]][,2]), col = "black", lwd = 2, lty = 3)
      lines(tseq2, as.numeric(pred[[1]][,3]), col = "black", lwd = 2, lty = 3)
    } else {
      points(tseq, as.numeric(alcohol[,2]), col = "red")
      lines(tseq1, as.numeric(trend[,2]), lwd = 2, col = "red")
      lines(tseq2, as.numeric(pred[[2]][,1]), col = "black", lwd = 3, lty = 1)
      lines(tseq2, as.numeric(pred[[2]][,2]), col = "black", lwd = 2, lty = 3)
      lines(tseq2, as.numeric(pred[[2]][,3]), col = "black", lwd = 2, lty = 3)
    }
    legend("topleft",
           legend = c(
             paste(colnames(alcohol)[1],sep="", " data"),
             paste(colnames(alcohol)[2],sep="", " data"),
             "filtered", "predictions", "95% C.I."
             ),
           pch = c(1,1,NA,NA,NA),
           lty = c(NA,NA,1,1,3),
           lwd = c(NA,NA,2,3,2),
           col = c("blue","red", "black", "black", "black"),
           bty = "n")
  }
  browser()


############################################################
  data("alcohol")
  alcoholPred          <- window(alcohol, start = 1969, end = 2007)
  par(mfrow = c(1,1), pty = "s")
  ts.plot(alcohol[,1:4], lty = c(1,2,3,4), ylab = "deaths")
  legend("topleft", lty= c(1,2,3,4), legend = colnames(alcohol)[1:4])
  model_alc_gauss4  <- SSModel(alcoholPred[,1:4] ~
                                 SSMtrend(2, Q = list(matrix(NA,4,4), matrix(0,4,4))) +
                                 SSMcustom(Z = diag(1,4), T = diag(0,4), Q = matrix(NA,4,4),
                                           P1 = matrix(NA,4,4)), distribution = "gaussian",
                               u = alcoholPred[,5:8]
  )
  updatefn <- function(pars, model, ...) {
    Q <- diag(pars[1:4])
    Q[upper.tri(Q)] <- pars[5:10]
    model["Q", etas = "level"] <- crossprod(Q)
    Q <- diag(pars[11:14])
    Q[upper.tri(Q)] <- pars[15:20]
    model["Q", etas = "custom"] <- model["P1", states = "custom"] <- crossprod(Q)
    model
  }
  init     <- chol(cov(alcoholPred[,1:4]))
  fitinit  <- fitSSM(model_alc_gauss4, updatefn = updatefn,
                     inits = rep(c(diag(init), init[upper.tri(init)]),2),
                     method = "BFGS")
  -fitinit$optim.out$val

  fit <- fitSSM(model_alc_gauss4, updatefn = updatefn,
                inits = fitinit$optim.out$par,
                method = "BFGS", nsim = 250)
  -fitinit$optim.out$val
  varcor <- fit$model["Q", etas = "level"]
  varcor[upper.tri(varcor)] <- cov2cor(varcor)[upper.tri(varcor)]
  print(varcor, digits = 2)

  varcor <- fit$model["Q", etas = "custom"]
  varcor[upper.tri(varcor)] <- cov2cor(varcor)[upper.tri(varcor)]
  print(varcor, digits = 2)
  out <- KFS(fit$model, nsim = 1000)
  print(out)
  plot(coef(out, states = c("level", "custom")), main = "Smoothed states",yax.flip = TRUE)
  browser()
  res <- rstandard(KFS(fit$model))
  browser()
  acf(res, na.action = na.pass)
  browser()
  pred <- predict(fit$model,
                  newdata <- SSModel(ts(matrix(NA,6,4), start = 2008) ~ -1 +
                                       SSMcustom(Z = fit$model$Z, T = fit$model$T, R = fit$model$R, Q = fit$model$Q),
                                     u = 1, distribution = "gaussian"),
                  interval = "confidence", nsim = 10000)
  trend <- signal(out, "trend")$signal
  par(mfrow = c(2,2), mar = c(2,2,2,2) + 0.1, oma = c(2,2,0,0))
  for(i in 1:4) {
    ts.plot(alcohol[,i], trend[,i], pred[[i]],
            xlab = NULL, ylab = "deaths",
            main = colnames(alcohol)[i])
  }
  browser()



###########################################################################
  data("alcohol")
  deaths     <- window(alcohol[,2], end = 2007)
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
  par(mfrow = c(1,1), pty = "s")
  plot(deaths/population, col = gray(0.5), typ = "b")
  model_gaussian <- SSModel(deaths/population ~ -1 + SSMcustom(Z = Zt, T = Tt,
              R = Rt, Q = Qt, a1 = a1, P1 = P1, P1inf = P1inf), H = Ht)
  fit_gaussian   <- fitSSM(model_gaussian, inits = c(0,0), method = "BFGS")
  out_gaussian   <- KFS(fit_gaussian$model)
  attributes(out_gaussian)
  lines(out_gaussian$att[,1], col = "red", lwd = 2)
  lines(out_gaussian$muhat, col = "blue", lwd = 2)
  title(main = "Alcohol", sub = "Gaussian")
  legend("topleft", legend = c("Observations", "One-step head pred.", "Smoothed estimates"),
  lty = c(NA,1,1), col = c(gray(0.5),"red","blue"), lwd = c(NA,2,2), pch = c(1,NA,NA), bty = "n"
  )
  browser()

############################################################
  data(Nile)
  level <- Nile
  par(mfrow = c(1,1), pty = "s")
  plot(level, typ = "b", col = gray(0.5))
  model_nile <- SSModel(level ~ -1 +
                  SSMcustom(Z = Zt, T = Tt,
                  R = Rt, Q = Qt, a1 = a1, P1 = P1, P1inf = P1inf), H = Ht)
  fit_nile   <- fitSSM(model_nile, inits = c(0,0), method = "BFGS")
  out_nile   <- KFS(fit_nile$model)
  attributes(out_nile)
  lines(out_nile$att[,1], col = "red", lwd = 2)
  lines(out_nile$muhat, col = "blue", lwd = 2)
  title(main = "Nile", sub = "Gaussian")
  legend("topleft", legend = c("Observations", "One-step head pred.", "Smoothed estimates"),
         lty = c(NA,1,1), col = c(gray(0.5),"red","blue"), lwd = c(NA,2,2), pch = c(1,NA,NA), bty = "n"
  )
#  pred_nile <- predict(model_nile,  interval = "prediction", level = 0.95)
  browser()


############################################################
  data("alcohol")
  alcoholPred        <- window(alcohol, start = 1969, end = 2007)
  ts.plot(alcohol[,1:4]/alcohol[,5:8], lty = c(1,2,3,4), ylab = "death rate")
  legend("topleft", lty= c(1,2,3,4), legend = colnames(alcohol)[1:4])
  model_alcohol  <- SSModel(alcoholPred[,1:4] ~
                        SSMtrend(2, Q = list(matrix(NA,4,4), matrix(0,4,4))) +
                        SSMcustom(Z = diag(1,4), T = diag(0,4), Q = matrix(NA,4,4),
                            P1 = matrix(NA,4,4)), distribution = "poisson",
                            u = alcoholPred[,5:8]
  )
  updatefn <- function(pars, model, ...) {
    Q <- diag(exp(pars[1:4]))
    Q[upper.tri(Q)] <- pars[5:10]
    model["Q", etas = "level"] <- crossprod(Q)
    Q <- diag(exp(pars[11:14]))
    Q[upper.tri(Q)] <- pars[15:20]
    model["Q", etas = "custom"] <- model["P1", states = "custom"] <- crossprod(Q)
    model
  }
  init     <- chol(cov(log(alcoholPred[,1:4]/alcoholPred[,5:8]))/10)
  fitinit  <- fitSSM(model_alcohol, updatefn = updatefn,
                     inits = rep(c(log(diag(init)), init[upper.tri(init)]),2),
                     method = "BFGS")
  -fitinit$optim.out$val

  fit <- fitSSM(model_alcohol, updatefn = updatefn,
                inits = fitinit$optim.out$par,
                method = "BFGS", nsim = 250)
  -fitinit$optim.out$val
  varcor <- fit$model["Q", etas = "level"]
  varcor[upper.tri(varcor)] <- cov2cor(varcor)[upper.tri(varcor)]
  print(varcor, digits = 2)

  varcor <- fit$model["Q", etas = "custom"]
  varcor[upper.tri(varcor)] <- cov2cor(varcor)[upper.tri(varcor)]
  print(varcor, digits = 2)
  out <- KFS(fit$model, nsim = 1000)
  print(out)
  browser()
  plot(coef(out, states = c("level", "custom")), main = "Smoothed sates",yax.flip = TRUE)
  browser()
  res <- rstandard(KFS(fit$model, filtering = "mean", smoothing = "none", nsim = 1000))
  acf(res, na.action = na.pass)
  pred <- predict(fit$model,
            newdata <- SSModel(ts(matrix(NA,6,4), start = 2008) ~ -1 +
            SSMcustom(Z = fit$model$Z, T = fit$model$T, R = fit$model$R, Q = fit$model$Q),
                  u = 1, distribution = "poisson"),
                  interval = "confidence", nsim = 10000)
  trend <- exp(signal(out, "trend")$signal)
  par(mfrow = c(2,2), mar = c(2,2,2,2) +0.1, oma = c(2,2,0,0))
  for(i in 1:4) {
    ts.plot(alcohol[,i]/alcohol[,4 + i], trend[,i], pred[[i]],
            xlab = NULL, ylab = "death rate",
            main = colnames(alcohol)[i])
  }
  browser()

###########################################################

}




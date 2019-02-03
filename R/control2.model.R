#' \code{control2.model} is a wrapper for controlling speed to obtain a target speed.
#'
#' @param tdf1df11 speed and location data, data.frame
# #' @examples
# #' control2.model()
#' @usage control2.model(tdf1df11)
#' @export
control2.model <- function(tdf1df11) {
  df     <- tdf1df2[,c(2,3,5,6)]
  start  <- 0
  end    <- 40
  tseq   <- seq(start, end, 0.125)
  u      <- round(53*5280/3600,0)
  usd    <- round(5280/3600*5,1)
  data   <- ts(df, start, end, frequency = 8)
  u0     <- 53*5280/3600
  u0     <- rep(u0, dim(df)[1])
  h      <- df[,2] - df[,4]
  htarget<- hsafe(u0[1],14)
  ### Differencing ##################################################################
  data   <- ts(data.frame(u0, u = df[,1],  times = tseq, u2 = df[,3], h, htarget),
               start, end, frequency = 8)
  par(mfrow = c(2,2), pty = "s")
  ts.plot(data[,1:2],
          xlab = "t, seconds",
          ylab = expression("Speed: "*dot(x)[t]),
          col = c("gold", "black"),
          lwd = c(6,2),
          ylim = c(40,140),
          main = "Lead Vehicle"
  )
  abline(h = 0,  col = gray(0.5))
  abline(v = 0,  col = gray(0.5))

  ts.plot(diff(data[,2]),
          xlab = "t, seconds",
          col = "black",
          ylim = c(-10,10),
          ylab = expression("Acceleration: "*nabla*dot(x)[t])
  )
  abline(h = 0, col = gray(0.5))
  abline(v = 0, col = gray(0.5))


  acf(diff(data[,2]), main = "")
  title(main = expression(nabla*dot(x)[t]), sub = "acf")

  pacf(diff(data[,2]),  main = "")
  title(main = expression(nabla*dot(x)[t]), sub = "pacf")

  browser()

### Lead Vehicle ###########################################################################################
  model <- SSModel(data[,2] ~ -1 +
                     SSMcustom(Z = matrix(c(1, 0), 1, 2),
                               T = array(diag(2), c(2, 2, nrow(data))),
                               Q = array(0, c(2, 2, nrow(data))),
                               P1inf = diag(2), P1 = diag(0, 2)), data = data)
  model$T[1, 2, ] <- c(diff(data[,3]), 1)
  model$Q[1, 1, ] <- c(diff(data[,3]), 1)^3/3
  model$Q[1, 2, ] <- model$Q[2, 1, ] <- c(diff(data[,3]), 1)^2/2
  model$Q[2, 2, ] <- c(diff(data[,3]), 1)

  updatefn <- function(pars, model, ...){
    model["H"] <- exp(pars[1])
    model["Q"] <- model["Q"] * exp(pars[2])
    model
  }
  fit <- fitSSM(model, inits = c(4, 4), updatefn = updatefn, method = "BFGS")

  pred <- predict(fit$model, interval = "prediction", level = 0.95)
  pred2 <- predict(fit$model, interval = "confidence", level = 0.95)
  par(mfrow = c(2,2), pty = "s")
  plot(x = data[,3], y = data[,2],
       pch = 16, cex = 0.75,
       main = "Lead Vehicle",
       col  = gray(0.6),
       ylim = c(40,140),
       xlab = "t, seconds",
       ylab = expression(dot(x)[t]))
  lines(x = tseq, y = u0,  lwd = 6, col = "gold")
  lines(x = tseq, y = pred[, 1], lwd = 3, col = "blue")
  lines(x = tseq, y = pred[, 2], lwd = 1, lty = 1, col = "blue")
  lines(x = tseq, y = pred[, 3], lwd = 1, lty = 1, col = "blue")
  abline(h = 0, col = gray(0.5))
  abline(v = 0, col = gray(0.5))

  legend("topleft",
         legend = c("Target", "Observed", "Predicted", "95% prediction interval" ),
         pch = c(NA,16,NA,NA),
         lty = c(1,NA,1,1),
         lwd = c(6,NA,3,1),
         col = c("gold", gray(0.5), "blue","blue"),
         bty = "n")

  ### Following Vehicle ###########################################################################################
  model <- SSModel(data[,4] ~ -1 +
                     SSMcustom(Z = matrix(c(1, 0), 1, 2),
                               T = array(diag(2), c(2, 2, nrow(data))),
                               Q = array(0, c(2, 2, nrow(data))),
                               P1inf = diag(2), P1 = diag(0, 2)), data = data)
  model$T[1, 2, ] <- c(diff(data[,3]), 1)
  model$Q[1, 1, ] <- c(diff(data[,3]), 1)^3/3
  model$Q[1, 2, ] <- model$Q[2, 1, ] <- c(diff(data[,3]), 1)^2/2
  model$Q[2, 2, ] <- c(diff(data[,3]), 1)
  updatefn <- function(pars, model, ...){
    model["H"] <- exp(pars[1])
    model["Q"] <- model["Q"] * exp(pars[2])
    model
  }
  fit <- fitSSM(model, inits = c(4, 4), updatefn = updatefn, method = "BFGS")

  pred <- predict(fit$model, interval = "prediction", level = 0.95)
  pred2 <- predict(fit$model, interval = "confidence", level = 0.95)

  plot(x = data[,3], y = data[,4],
       pch = 16, cex = 0.75,
       main = "Following Vehicle",
       col  = gray(0.6),
       ylim = c(40,140),
       xlab = "t, seconds",
       ylab = expression(dot(x)[t]))
  lines(x = tseq, y = u0,  lwd = 6, col = "gold")
  lines(x = tseq, y = pred[, 1], lwd = 3, col = "red")
  lines(x = tseq, y = pred[, 2], lwd = 1, lty = 1, col = "red")
  lines(x = tseq, y = pred[, 3], lwd = 1, lty = 1, col = "red")
  abline(h = 0, col = gray(0.5))
  abline(v = 0, col = gray(0.5))

  legend("topleft",
         legend = c("Target", "Observed", "Predicted", "95% prediction interval" ),
         pch = c(NA,16,NA,NA),
         lty = c(1,NA,1,1),
         lwd = c(6,NA,3,1),
         col = c("gold", gray(0.5), "red","red"),
         bty = "n")
  browser()
  ### Headway ##########################################################################################

  plot(tseq,tdf1df2[,3], lwd = 2, typ = "l",
       col = "blue",
       xlab = "t, seconds",
       ylab = expression(x[t]))
  lines(tseq,tdf1df2[,6], lwd = 2, col = "red")
  abline(h = c(0,-500), col = gray(0.5))
  abline(v = 0, col = gray(0.5))
  title(main = "Bottleneck Merge")
  legend("topleft",
         legend = c("Lead vehicle", "Following vehicle" ),
         lty = c(1,1),
         lwd = c(1,1),
         col = c("blue","red"),
         bty = "n")

  model <- SSModel(data[,5] ~ -1 +
                     SSMcustom(Z = matrix(c(1, 0), 1, 2),
                               T = array(diag(2), c(2, 2, nrow(data))),
                               Q = array(0, c(2, 2, nrow(data))),
                               P1inf = diag(2), P1 = diag(0, 2)), data = data)
  model$T[1, 2, ] <- c(diff(data[,3]), 1)
  model$Q[1, 1, ] <- c(diff(data[,3]), 1)^3/3
  model$Q[1, 2, ] <- model$Q[2, 1, ] <- c(diff(data[,3]), 1)^2/2
  model$Q[2, 2, ] <- c(diff(data[,3]), 1)

  updatefn <- function(pars, model, ...){
    model["H"] <- exp(pars[1])
    model["Q"] <- model["Q"] * exp(pars[2])
    model
  }
  fit <- fitSSM(model, inits = c(4, 4), updatefn = updatefn, method = "BFGS")

  pred <- predict(fit$model, interval = "prediction", level = 0.95)
  pred2 <- predict(fit$model, interval = "confidence", level = 0.95)
  plot(x = as.numeric(data[,3]),
          y = as.numeric(data[,5]),
       pch = 1, cex = 1.5,
       main = "Headways",
       col  = gray(0.6),
       xlab = "t, seconds",
       ylab = expression(h[t]))
  lines(x = tseq, y = data[,6],  lwd = 6, col = "gold")
  lines(x = tseq, y = pred[, 1], lwd = 3, col = "red")
  lines(x = tseq, y = pred[, 2], lwd = 1, lty = 1, col = "red")
  lines(x = tseq, y = pred[, 3], lwd = 1, lty = 1, col = "red")
  abline(h = 0, col = gray(0.5))
  abline(v = 0, col = gray(0.5))

  legend("topleft",
         legend = c("Target", "Observed", "Predicted", "95% prediction interval" ),
         pch = c(NA,1,NA,NA),
         lty = c(1,NA,1,1),
         lwd = c(6,NA,3,1),
         col = c("gold", gray(0.5), "red","red"),
         bty = "n")
  browser()
}

#' \code{control1.model} is a wrapper for controlling speed to obtain a target speed.
#'
# #' @examples
# #' control1.model()
#' @usage control1.model()
#' @export
control1.model <- function() {
  set.seed(123)
  start  <- 0
  end    <- 40
  tseq   <- seq(start, end, 0.125)
  u0     <- round(20*5280/3600,0)
  a      <- 2
  u      <- a * tseq
  df1    <- data.frame(u = rep(20,length(u[u <= 20])))
  df2    <- data.frame(u = u[u > 20])
  df     <- rbind(df1,df2)
  usd    <- rnorm(dim(df)[1],0,5*5280/3600)
  u      <- df$u + usd
  u0     <- rep(40, dim(df)[1])
  par(mfrow = c(1,2), pty = "s")
  acf(diff(u), main = "")
  title(main = expression("Arima, "*nabla*dot(x)[t]), sub = "acf")
  pacf(diff(u), main = "")
  title(main = expression("Arima, "*nabla*dot(x)[t]), sub = "pacf")
  browser()

  #####################################################################
  data   <- ts(data.frame(u0, u, times = tseq), start, end, frequency = 8)
  par(mfrow = c(1,1), pty = "s")
  ts.plot(data[,1:2], xlab = "t, seconds", ylab = expression("Speed: "*dot(x)[t]),
          col = c("gold", "black"),
          lwd = c(6,2),
          ylim = c(-10,100)
  )
  abline(h = 0,  col = gray(0.5))
  abline(v = 0,  col = gray(0.5))
  browser()

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

  plot(x = data[,3], y = data[,2], pch = 19, col = gray(0.6),
       xlab = "t, seconds", ylab = "u(t), feet per second")
  lines(x = tseq, y = pred[, 1], lwd = 3, col = "blue")
  lines(x = tseq, y = pred[, 2], lwd = 3, lty = 2, col = "blue")
  lines(x = tseq, y = pred[, 3], lwd = 3, lty = 2, col = "blue")
  lines(x = tseq, y = pred2[, 2], lwd = 2, lty = 2, col = "orange")
  lines(x = tseq, y = pred2[, 3], lwd = 2, lty = 2, col = "orange")

 browser()


  ### model ############################################################
  model  <- SSModel(data[,2] ~  -1 +
                      SSMcustom(Z  = matrix(c(0,1,0),1,3),
                                T  = matrix(c(1,0,0,0.125,1,0,0.125^2/2,0.125,1),3,3),
                                R  = matrix(c(0,0,0,0,1,0,0,0,0),3,3),
                                Q  = matrix(c(NA,0,0,0,NA,0,0,0,NA),3,3),
                                a1 = matrix(c(0,1,0),1,3),
                                P1 = matrix(0,3,3),
                                P1inf = diag(3)
                      ),
                    distribution = "gaussian",
                    tol  = .Machine$double.eps^0.5,
                    H = NA
  )

  ###########################################################################################################
  check_model  <- function(model) (    model$H[1,1,1] > 0 &
                                       model$Q[1,1,1] > 0 &
                                       model$Q[2,2,1] > 0 &
                                       model$Q[3,3,1] > 0
  )
  update_model <- function(pars, model) {
    model["Q"][1,1,1] <- model["Q"][1,1,1]*exp(pars[1])
    model["Q"][2,2,1] <- model["Q"][2,2,1]*exp(pars[2])
    model["Q"][3,3,1] <- model["Q"][3,3,1]*exp(pars[3])
    model["H"]              <- exp(pars[4])
    model
  }
  browser()

  fit     <- fitSSM(model,
                    inits    = c(4,4,4,4),
                    updatefn = update_model,
                    checkfn  = check_model,
                    method   = "BFGS")
  out     <- KFS(fit$model)

  Out <- ts(data.frame(stand = rep(u,321),
                       obs = data[,1],
                       smooth = signal(out)$signal[[1]],
                       prd    = signal(out,filtered = TRUE)$signal),
            start, end, frequency = 8)
  par(mfrow = c(1,1), pty = "s")
  ts.plot(window(Out, start = c(0,1), end = c(40,0)),
          col = c("gold", gray(0.5), "black", "blue"),
          ylim = c(0,150),
          ylab = "u(t), feet per second", lty = c(1,3,1,1), lwd = c(6,2,2,2)
  )
  title(main = expression(dot(x)[t]))
  legend("bottomright",
         legend = c("target", "Observed", "One-step ahead predictions", "Smoothed estimates" ),
         lty = c(1,3,1,1),
         lwd = c(6,2,2,2),
         col = c("gold", gray(0.5), "blue","black"),
         bty = "n")

  legend("topleft",c(
    expression(""),
    bquote(bar(u) == .(u)),
    bquote(sigma[w] == .(usd))),
    bty = "n"
  )



  browser()

}

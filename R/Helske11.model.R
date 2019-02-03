#' \code{Helske2.model} is a wrapper for estimating speed and location with Kalman filtering \code{SSMcustom} See \code{Helske9.R}.
#'
#' @param tdf1df11 speed and location data, data.frame
# #' @examples
# #' Helske11.model(tdf1df2)
#' @usage Helske11.model(tdf1df2)
#' @export
Helske11.model <- function(tdf1df2) {
  df     <- tdf1df2[,c(2,3,5,6)]
  start  <- 0
  end    <- 40
  u      <- round(53*5280/3600,0)
  usd    <- round(5280/3600*5,1)
  data   <- ts(df, start, end, frequency = 8)
  par(mfrow = c(1,2), pty = "s")
  acf(diff(df[,1]))
  title(main = expression("Arima, "*nabla*dot(x)[t]), sub = "acf")

  pacf(diff(df[,1]),  main = expression(nabla*dot(x)[t]))
  title(main = expression("Arima, "*nabla*dot(x)[t]), sub = "pacf")
  browser()
  y  <- ts(diff(df[,1]), start = c(0,1), end = 40, frequency = 8)
  ts.plot(y,
          col = "black",
          ylim = c(-10,10)
  )

  browser()
#############################################################################################################
  u      <- mean(data[,1])
  par(mfrow = c(1,1), pty = "s")

  drift <- 1:length(y)
  model_arima <- SSModel(data[,1] ~ drift +
                           SSMarima(ar = c(0, 0), d = 1, ma = c(0, 0), Q = 1) +
                           SSMcustom(Z  = matrix(c(1,0.125),1,2),
                                     T  = matrix(c(1,0,1,1),2,2),
                                     Q  = matrix(NA,1,1),
                                     P1 = matrix(c(1,0,0,1),2,2)
                           ),
                         distribution = "gaussian",
                         tol  = .Machine$double.eps^0.5
  )
  print(model_arima)
  update_model <- function(pars, model) {
    tmp <- SSMarima(ma = artransform(pars[1:2]),
                    d = 1,
                    ar = artransform(pars[3:4]),
                    Q = exp(pars[5]))
    model["T", states = "arima"] <- tmp$T
    model["R", states = "arima"] <- tmp$R
    model["Q", states = "arima"] <- tmp$Q
    model["P1",states = "arima"] <- tmp$P1
    model["Q", states = "custom"] <- exp(pars[6])
    model["P1",states = "custom"] <- exp(pars[7])
    model
  }
  fit_arima <- fitSSM(model_arima, inits = rep(1,7), updatefn = update_model,
                      method = "L-BFGS-B", lower = c(-1, 0), upper = c(1, 100))
  print(fit_arima$optim.out$par)
  out_arima <- KFS(fit_arima$model)
  Out   <- ts(data.frame(u, obs = data[,1],
                         smooth = signal(out_arima)$signal[[1]],
                         prd    = signal(out_arima,filtered = TRUE)$signal
  ), start = start, end = end, frequency = 8)


  ###########################################################################################################

  par(mfrow = c(1,1), pty = "s")
  ts.plot(window(Out, start = c(0,1), end = c(40,0)),
          col = c("gold", gray(0.5), "black", "blue"),
          ylim = c(0,150),
          ylab = "u(t), feet per second", lty = c(1,3,1,1), lwd = c(6,2,2,2)
  )
  title(main = expression(dot(x)[t]))
  legend("bottomright",
         legend = c("Target", "Observed", "Filtered", "Signal" ),
         lty = c(1,3,1,1),
         lwd = c(6,2,2,2),
         col = c("gold", gray(0.5), "blue","black"),
         bty = "n")
  u = round(u)
  legend("topleft",c(
    expression(""),
    bquote(bar(u) == .(u)),
    bquote(sigma[U] == .(usd))),
    bty = "n"
  )

  browser()

}

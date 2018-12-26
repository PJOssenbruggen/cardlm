#' \code{xuForecast} Random walk model for forecasting.
#'
#' All observations are in units of feet and seconds.
#' @param tdf1fd2, observations from \code{cartools}, a matrix.
#' @param veh, vehicle number.
#' @param dV, variance of the observation noise.
#' @param dW, diagonal elements of the variance matrix of the system noise.
#' @param m0, the expected value of the pre-sample state vector.
#' @param C0, the variance matrix of the pre-sample state vector.
#' @param tend, end time in seconds
#' @param frequency, number of time-steps per second
#' @param wv signal noise ratio, a number
#' @return A \code{ts} plot to study the effects of \code{dW} on the forecast.
#' @usage xuForecast(tdf1df2, veh, dV, dW, m0, C0, tend, frequency, wv)
# #' @examples
# #' xuForecast(tdf1df2, 1, 7.33, 7.33, 77.88, 77.88, 40, 8, 0.01)
#' @export
xuForecast <- function(tdf1df2, veh, dV, dW, m0, C0, tend, frequency, wv) {
  N       <- as.numeric(dim(tdf1df2)[1])
  u       <- as.matrix(tdf1df2[,seq(2,dim(tdf1df2)[2],3)])
  x       <- as.matrix(tdf1df2[,seq(3,dim(tdf1df2)[2],3)])
  u       <- u[,veh]
  x       <- x[,veh]
  xu      <- cbind(x,u)
  y       <- ts(xu, start = 0, end = tend,  frequency)
  tseq    <- seq(0,tend,1/frequency)
  lst     <- xuMLEparam(tdf1fd2, veh, tend, frequency)
  FF      <- lst[[3]]
  m0      <- as.numeric(y[1,])
  C0      <- lst[[2]]
  GG      <- lst[[4]]
  V       <- lst[[5]]
  W       <- lst[[6]]
  W[1,2]  <- W[2,1] <- 0
  xu      <- dlm::dlm(m0 = m0, C0 = C0, FF = FF, GG = GG, V = V, W = W)
  Foremod <- dlm::dlmFilter(y, xu)
  a       <- as.data.frame(cbind(t = tseq, y, Foremod$f))
  #  str(Foremod,1)
  par(mfrow = c(1,2), pty = "s")
  plot(a[,1],a[,3],xlab = "t, seconds", ylab = "u, feet per second",
        ylim = c(0,1.1*max(a[,3])), typ = "l", lwd = 4, col = gray(0.0))
  lines(a[,1], a[,5], lty = 1, col = "orange", lwd = 1)
  title(main = "Kalman Forecasting")
  axis(side = 4, at = m0[2], expression(bar(u)))
  abline(h = m0[2], col = gray(0.8))
  abline(v = 0, col = gray(0.8))
  wvratio <- W[2,2]/V[2,2]
  legend("topleft", legend = c(paste("data, V = ", round(V[2,2],4)),
                               paste("one-step ahead, W = ", round(W[2,2],2))),
         lty = c(1,1), col = c(gray(0.5), "orange"),
         lwd = c(3,1), bty = "n")
  plot(a[,1],a[,2],xlab = "t, seconds", ylab = "x, feet", typ = "l",lwd = 4)
  lines(a[,1], a[,4], lty = 1, col = "orange", lwd = 1)
  title(main = "Kalman Forecasting")
  abline(h = c(0), col = gray(0.8))
  abline(v = 0, col = gray(0.8))
  wvratio <- W[1,1]/V[1,1]
  legend("topleft", legend = c(paste("data, V = ", round(V[1,1],4)),
                               paste("one-step ahead, W = ", round(W[1,1],2))),
         lty = c(1,1), col = c(gray(0.5), "orange"),
         lwd = c(3,1), bty = "n")
}

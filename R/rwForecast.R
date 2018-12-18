#' \code{rwForecast} Random walk model for forecasting.
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
#' @usage rwForecast(tdf1df2, veh, dV, dW, m0, C0, tend, freqency, wv)
#' @examples
#' rwForecast(tdf1df2, 1, 7.33, 7.33, 77.88, 77.88, 40, 8, 0.01)
#' @export
rwForecast <- function(tdf1df2, veh, dV, dW, m0, C0, tend, freqency, wv) {
  N     <- as.numeric(dim(tdf1df2)[1])
  u     <- as.matrix(tdf1df2[,seq(2,dim(tdf1df2)[2],3)])
  u     <- u[,veh]
  u.ts  <- ts(u, start = 0, end = tend,  frequency)
  plot(u.ts, xlab = "t, seconds", ylab = "u, feet per second",
       ylim = c(0,1.1*max(tdf1df2[,2])))
  rw        <- dlm::dlmModPoly(order = 1, dV, dW, C0, m0)
  wvratio   <- round(dW/dV, 2)
  print(unlist(rw))
  Foremod   <- dlm::dlmFilter(u, rw)
#  str(Foremod,1)
  a         <- cbind(u.ts,Foremod$f)
  lines(a[,2], lty = 2, col = "blue", lwd = 2)
  title(main = "Kalman Forecasting")
  axis(side = 4, at = m0, expression(bar(u)))
  abline(h = m0, col = gray(0.8))
  dW      <- wv * dW
  wvratio2<- round(dW/dV, 6)
  rw2     <- dlm::dlmModPoly(order = 1, dV, dW, C0, m0)
  print(unlist(rw2))
  Foremod2 <- dlm::dlmFilter(u, rw2)
#  str(Foremod2,1)
  a         <- cbind(u.ts,Foremod$f, Foremod2$f)
  lines(a[,3], lty = 1, col = "red", lwd = 2)
  legend("topleft", legend = c(paste("data, V = ", dV),
                               paste("one-step ahead, W/V = ", wvratio),
                               paste("one-step ahead, W/V = ", wvratio2)),
         lty = c(1,3), col = c(gray(0.5), "blue"),
         lwd = c(2,2), bty = "n")
}

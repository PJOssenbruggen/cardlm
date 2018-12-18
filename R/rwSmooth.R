#' \code{rwSmooth} Random walk model for smoothing.
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
#' @usage rwSmooth(tdf1df2, veh, dV, dW, m0, C0, tend, freqency)
#' @examples
#' rwSmooth(tdf1df2, 1, 73.3, 73.3, 77.88, 77.88, 40, 8)
#' @export
rwSmooth <- function(tdf1df2, veh, dV, dW, m0, C0, tend, freqency) {
  N     <- as.numeric(dim(tdf1df2)[1])
  u     <- as.matrix(tdf1df2[,seq(2,dim(tdf1df2)[2],3)])
  u     <- u[,veh]
  u.ts  <- ts(u, start = 0, end = tend,  frequency)
  plot(u.ts, xlab = "t, seconds", ylab = "u, feet per second",
       ylim = c(0,1.1*max(tdf1df2[,2])))
  rw        <- dlm::dlmModPoly(order = 1, dV, dW, C0, m0)
  wvratio   <- round(dW/dV, 2)
  print(unlist(rw))
  title(main = "Kalman Smoothing")
  axis(side = 4, at = m0, expression(bar(u)))
  abline(h = m0, col = gray(0.8))
  Smthmod   <- dlm::dlmSmooth(u.ts, rw)
  attach(Smthmod)
  a         <- cbind(u.ts, Smthmod$s)
  lines(a[,2], lty = 2, col = "orange", lwd = 2)
  N         <- as.numeric(dim(tdf1df2)[1])
  class(U.S[[N+1]])
  drop(dlm::dlmSvd2var(U.S[[N+1]], D.S[[N+1]]))
  unlist(dlm::dlmSvd2var(U.S,D.S))
  hwid     <- qnorm(0.025, lower.tail = FALSE) * sqrt(unlist(dlm::dlmSvd2var(U.S,D.S)))
  detach(Smthmod)
  a        <- cbind(a, as.vector(Smthmod$s) + hwid %o% c(-1,1))
  lines(a[,3], lty = 1, col = "orange", lwd = 1)
  lines(a[,4], lty = 1, col = "orange", lwd = 1)
  legend("topleft", legend = c(paste("data, V = ", dV),
                               paste("smoothed, W/V = ", wvratio),
                               "95% C.I."),
         lty = c(1,1,1), col = c(gray(0.5), "orange"),
         lwd = c(1,2,1), bty = "n")
  return(u)
}

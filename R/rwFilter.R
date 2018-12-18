#' \code{rwFilter} Random walk model for filtering and smoothing.
#'
#' @param tdf1fd2, \code{cartools} data set, a matrix.
#' @return Plots of the speed \code{u} data, Kalman filter and smoother.
#' @usage rwFilter(tdf1df2)
#' @export
rwFilter <- function(tdf1df2) {
  plot(tdf1df2[,1], tdf1df2[,2], typ = "b", col = gray(0.5), xlab = "t, seconds", ylab = "u, feet per second",
       ylim = c(0,1.1*max(tdf1df2[,2])))
  rw  <- dlm::dlmModPoly(order = 1, dV = 73.3, dW = 7.73, C0 = 77.88, m0 = 77.88)
  snr <- 0.1
  print(unlist(rw))
  ulead     <- tdf1df2[,2]
  Filtmod   <- dlm::dlmFilter(ulead, rw)
  lines(tdf1df2[,1], dlm::dropFirst(Filtmod$m), lwd = 2, col = "blue")
  Smthmod <- dlm::dlmSmooth(ulead, rw)
  N <- as.numeric(dim(tdf1df2)[1])
  str(Smthmod,1)
  print(N)
  attach(Smthmod)
  class(U.S[[N+1]])
  drop(dlm::dlmSvd2var(U.S[[N+1]], D.S[[N+1]]))
  lines(tdf1df2[,1], dlm::dropFirst(Smthmod$s), lwd = 2, col = "orange")
  unlist(dlm::dlmSvd2var(U.S,D.S))
  hwid     <- qnorm(0.025, lower.tail = FALSE) * sqrt(unlist(dlm::dlmSvd2var(U.S,D.S)))
  smooth   <- cbind(s, as.vector(s) + hwid %o% c(-1,1))
  lines(tdf1df2[,1], dlm::dropFirst(Smthmod$s), lwd = 2, col = "orange")
  lines(tdf1df2[,1], dlm::dropFirst(Smthmod$s) + hwid[-1], lwd = 1, col = "orange")
  lines(tdf1df2[,1], dlm::dropFirst(Smthmod$s) - hwid[-1], lwd = 1, col = "orange")
  title(main = "Kalman Filtering")
  legend("topleft",legend = c("data", "filtered", "smoothed", "95% C.I."),
         lty = c(1,1,1,1), pch = c(1,NA,NA,NA), col = c(gray(0.5), "blue", "orange","orange"),
         lwd = c(2,2,2,1), bty = "n")
  abline(h=77.88, col = gray(0.8))
  axis(side = 4, at = 77.88, expression(bar(u)))
  axis(side = 3, at = tdf1df2[N/2,1], bquote(W/V == .(snr)), line = -1, tick = FALSE)
  detach("Smthmod")
}

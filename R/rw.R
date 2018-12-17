#' \code{rw} Random walk model.
#'
#' @param tdf1fd2, \code{cartools} data set, a matrix.
#' @return Plot of the data.
#' @usage rw(tdf1df2)
# #' @ examples
# #' rw(tdf1df2)
#' @export
rw <- function(tdf1df2) {
  plot(tdf1df2[,1], tdf1df2[,3], typ = "b", col = gray(0.5), xlab = "t, seconds", ylab = "x, feet")
  rw <- dlmModPoly(order = 1, dV = 500, dW = 5, C0 = -700, m0 = -700)
  print(unlist(rw))
  xlead    <- tdf1df2[,3]
  zip.rw.f <- dlmFilter(xlead, rw)
  zip.rw.s <- dlmSmooth(xlead, rw)
  attach(zip.rw.s)
  hwid     <- qnorm(0.025, lower.tail = FALSE) * sqrt(unlist(dlmSvd2var(U.S,D.S)))
  smooth   <- cbind(s, as.vector(s) + hwid %o% c(-1,1))
  lines(tdf1df2[,1], zip.rw.f$m[-1], lwd = 2, col = "blue")
  lines(tdf1df2[,1], zip.rw.s$s[-1], lwd = 2, col = "orange")
  lines(tdf1df2[,1], smooth[-1,3], lwd = 2, col = "tan")
  lines(tdf1df2[,1], smooth[-1,2], lwd = 2, col = "tan")
  detach(zip.rw.s)

}

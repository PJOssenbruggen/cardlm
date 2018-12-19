#' \code{MLEparam} MLE.
#'
#' @param tdf1fd2, observations from \code{cartools}, a matrix.
#' @param veh, vehicle number.
#' @param tend, end time in seconds
#' @param m0, the expected value of the pre-sample state vector.
#' @param C0, the variance matrix of the pre-sample state vector.
#' @param frequency, number of time-steps per second
#' @usage MLEparam(tdf1fd2, veh, tend, frequency, m0, C0)
# #' @examples
# #' MLEparam(tdf1fd2, 1, 40, 8, 77.88, 77.88)
#' @export
MLEparam <- function(tdf1fd2, veh, tend, frequency, m0, C0) {
  build <- function(param) {
    dlm::dlmModPoly(order = 1, dV = exp(parm[1]), dW = exp(parm[2]), m0, C0)
  }
  N     <- as.numeric(dim(tdf1df2)[1])
  u     <- as.matrix(tdf1df2[,seq(2,dim(tdf1df2)[2],3)])
  u     <- u[,veh]
  y     <- ts(u, start = 0, end = tend,  frequency)
  results <- {}
  for(i in 1:10) {
    parm  <- rep(i,2)
    fit   <- dlm::dlmMLE(y, parm, build, debug  = TRUE)
    Vest  <- unlist(build(fit$par))[c("V", "W")][1]
    West  <- unlist(build(fit$par))[c("V", "W")][2]
    results <- rbind(results, data.frame(i, V = Vest, W = West, value = fit$value))
  }
  print(results)
}

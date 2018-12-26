#' \code{xuMLEparam} MLE for \code{cartools} speed.
#'
#' @param tdf1fd2, observations from \code{cartools}, a matrix.
#' @param veh, vehicle number.
#' @param tend, end time in seconds
#' @param frequency, observation per second, a number
#' @usage xuMLEparam(tdf1fd2, veh, tend, frequency)
# #' @examples
# #' xuMLEparam(tdf1fd2, 1, 40, 8)
#' @export
xuMLEparam <- function(tdf1fd2, veh, tend, frequency) {
  N       <- as.numeric(dim(tdf1df2)[1])
  u       <- as.matrix(tdf1df2[,seq(2,dim(tdf1df2)[2],3)])
  x       <- as.matrix(tdf1df2[,seq(3,dim(tdf1df2)[2],3)])
  u       <- u[,veh]
  x       <- x[,veh]
  xu      <- cbind(x,u)
  y       <- ts(xu, start = 0, end = tend,  frequency)
  fit     <- dlm::dlmMLE(y, rep(0,6), xuBuild)
  lst     <- xuBuild(fit$par)
  return(lst)
}

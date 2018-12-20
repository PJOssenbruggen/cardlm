#' \code{rwMLEparam} MLE for \code{cartools} speed.
#'
#' @param tdf1fd2, observations from \code{cartools}, a matrix.
#' @param veh, vehicle number.
#' @param tend, end time in seconds
#' @param frequency, observation per second, a number
#' @usage rwMLEparam(tdf1fd2, veh, tend, frequency)
# #' @examples
# #' rwMLEparam(tdf1fd2, 1, 40, 8)
#' @export
rwMLEparam <- function(tdf1fd2, veh, tend, frequency) {
  N     <- as.numeric(dim(tdf1df2)[1])
  u     <- as.matrix(tdf1df2[,seq(2,dim(tdf1df2)[2],3)])
  u     <- u[,veh]
  y     <- ts(u, start = 0, end = tend,  frequency)
  parm  <- rep(0,2)
  fit   <- dlm::dlmMLE(y, rep(0,2), rwBuild)
  results <- unlist(rwBuild(fit$par))
  return(results)
}

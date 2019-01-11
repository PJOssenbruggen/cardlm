#' \code{Helske2.model} is a wrapper function for estimating a successful, safe merge.
#'
#' @usage uses \code{Helske2.model}, \code{xab} and \code{uab} from the \code{cartools} Package.
#' @param tdf1df2, time, speed and locations, a matrix
#' @param veh, a number
# #' @examples
# #' Helske2.model(tdf1df2, veh)
#' @export
Helske2.model <- function(tdf1df2, veh) {
  vehseq <- seq(2,31,3)
  z      <- tdf1df2[1:100, vehseq[veh]]
  df     <- data.frame(z)
  start  <- 0
  end    <- 100
  data   <- ts(df, start, end)
  model1 <- SSModel(z ~ SSMtrend(1, Q = NA), H = NA)
  fit1   <- fitSSM(model1, inits = c(0,0), method = "BFGS")
  out1   <- KFS(fit1$model)
  Qhat   <- round(as.numeric(fit1$model$Q),2)
  Hhat   <- round(as.numeric(fit1$model$H),2)
  df     <- data.frame(Z = as.numeric(fit1$model$Z),  Hhat,  Qhat)
  par(mfrow = c(1,2), pty = "s")
  ts.plot(ts(z), out1$a, out1$muhat,
          col = c(gray(0.5), "blue",  "black"),
          ylim = c(0.75 * min(ts(z)),1.1 * max(ts(z))),
          ylab = "u", lty = c(3,1,1), lwd = c(2,2,2)
  )
  title(main = "Predictions")
  legend("bottomright",
         legend = c(  "Observed", "One-step ahead", "Smoothed" ),
         lty = c(3,1,1),
         lwd = c(2,2,2),
         col = c(gray(0.5), "blue","black"),
         bty = "n")
  #  print(data.frame(ts(x), out1$a[-1], out1$muhat, ts(sim[,,1])))

  ts.plot(out1$P[1,,], col = "black", ylab = "P", lwd = 2)
  title(main = "Covariance")
  legend("topright",c(
    expression(""),
    bquote(hat(H) == .(Hhat)),
    bquote(hat(Q) == .(Qhat))),
    bty = "n"
  )

  ### Notes
  # 1. u = lead vehicle speed
  # 2. H = 0 suggests no noise in observations

  return(list(model1, fit1, out1, df, out1$v))
}

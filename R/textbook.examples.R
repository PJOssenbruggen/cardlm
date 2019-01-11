#' \code{textbook.examples} is a wrapper function for estimating a successful, safe merge.
#'
#' @usage uses \code{textbook.examples}, \code{xab} and \code{uab} from the \code{cartools} Package.
# #' @param tdf1df2, time, speed and locations, a matrix
# #' @examples
# #' textbook.examples()
#' @export
textbook.examples <- function() {
  ### Example 5.2 Maybeck page 220
  set.seed(1)
  start <- 0
  end   <- 10
  x     <- rep(1,10)
  # Added x to eliminate tsp error
  z     <- cumsum(rnorm(end, 0, sqrt(0.5)))
  df    <- data.frame(x, z)
  data  <- ts(df, start, end)
  Zt    <- matrix(1)
  Ht    <- matrix(NA)
  Tt    <- matrix(1)
  Rt    <- matrix(1)
  Qt    <- matrix(NA)
  a1    <- matrix(1)
  P1    <- matrix(0)
  P1inf <- diag(1)
  par(mfrow = c(1,1), pty = "s")
  model1 <- SSModel(z ~ -1 + SSMcustom(Z = Zt, T = Tt, R = Rt, Qt), H = Ht, data = data)
  fit1   <- fitSSM(model1, inits = c(0,0), method = "BFGS")
  out1   <- KFS(fit1$model)
  ts.plot(data[,1], data[,2], out1$a, out1$muhat,
          col = c("gold", gray(0.5), "blue",  "black"),
          ylab = "u(t), feet per second", lty = c(1,3,1,1), lwd = c(6,2,2,2)
  )
  abline(v = 0, col = gray(0.5))
  abline(h = 0, col = gray(0.5))
  title(main = "Gyro", sub = "Gaussian")
  legend("topleft",
         legend = c(expression(data), "Filtered estimates", "Smoothed estimates"),
         pch = c(1,NA,NA,NA),
         lty = c(1,1,1),
         lwd = c(1,4,2,3),
         col = c("blue","red", "yellow"),
         bty = "n")
}

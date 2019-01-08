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
  data  <- data.frame(z = cumsum(rnorm(end, 0, sqrt(0.5))))
  data  <- ts(data, start, end)
  Zt <- matrix(-1)
  Ht <- matrix(NA)
  Tt <- matrix(1)
  Rt <- matrix(1)
  Qt <- matrix(NA)
  a1 <- matrix(1)
  P1 <- matrix(0)
  P1inf <- diag(1)

  par(mfrow = c(1,1), pty = "s")
  model1 <- SSModel(z ~ -1 +
                SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, a1 = a1, P1 = P1, P1inf = P1inf),
                data = data, distribution = "poisson", u = data)
  print(model1$y)
  print(model1$Z)
  print(model1$H)
  print(model1$T)
  print(model1$R)
  print(model1$Q)
  browser()
  fit1   <- fitSSM(model1, inits = c(0,0), method = "BFGS")
  out1   <- KFS(fit1$model)
  browser()
  print(out1)
  print(out1$att)
  print(out1$a)
  print(out1$muhat)
  print(out1$Ptt)
  print(out1$P)
  par(mfrow = c(1,2), pty = "s")
  plot(data,type = "b", ylab = "y", xlab = "t, seconds",lwd = 2,col = "blue",ylim=c(-3,5) )
  lines(out1$att, col = "red", lwd = 4)
  lines(out1$muhat, col = "yellow", lwd = 2)
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
  plot(out1$Ptt,type = "l", ylab = "P", xlab = "t, seconds",lwd = 2,col = "black", ylim=c(0,1) )
  title(main = "Error Variance", sub = "Gaussian")
  abline(v = 0, col = gray(0.5))
  abline(h = 0, col = gray(0.5))
  legend("topleft",
         legend = c("a", "Ptt", "P"),
         pch = c(NA,NA,NA),
         lty = c(1,1,1),
         lwd = c(1,4,2),
         col = c("black","red", "yellow"),
         bty = "n")
  browser()
}

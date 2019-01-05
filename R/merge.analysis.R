#' \code{merge.analysis} is a wrapper function for estimating a successful, safe merge.
#'
#' @usage uses \code{merge.analysis}, \code{xab} and \code{uab} from the \code{cartools} Package.
#' @param tdf1df2, time, speed and locations, a matrix
# #' @examples
# #' merge.analysis()
#' @export
merge.analysis <- function(tdf1df2) {
  par(mfrow = c(1,1), pty = "s")
  plot(tdf1df2[,1], tdf1df2[,3], typ = "l", ylim = c(-1200,2400),
       ylab = expression("x, feet"), xlab = "t, seconds")
  xveh  <- seq(3,30,3)
  uveh  <- seq(2,29,3)
  tveh  <- tdf1df2[,1]
  h     <- x.xe <- u.x0  <- u.xe <- t.x0 <- t.xe <- rep(NA,10)
  u     <- tdf1df2[tdf1df2[,xveh[1]] <= 0, uveh[1]]
  t     <- tveh[tdf1df2[,xveh[1]] <= 0]
  u.x0[1] <- u[length(t)]
  t.x0[1] <- t[length(t)]
  x     <- tdf1df2[tdf1df2[,xveh[1]] <= -500, xveh[1]]
  u     <- tdf1df2[tdf1df2[,xveh[1]] <= -500, uveh[1]]
  t     <- tveh[tdf1df2[,xveh[1]] <= -500]
  x.xe[1] <- x[length(t)]
  u.xe[1] <- u[length(t)]
  t.xe[1] <- t[length(t)]
  for(i in 2:10) {
    lines(tdf1df2[,1], tdf1df2[,xveh[i]])
    u     <- tdf1df2[tdf1df2[,xveh[i]] <= 0, uveh[i]]
    t     <- tveh[tdf1df2[,xveh[i]] <= 0]
    u.x0[i] <- u[length(t)]
    t.x0[i] <- t[length(t)]
    x     <- tdf1df2[tdf1df2[,xveh[i]] <= -500, xveh[i]]
    u     <- tdf1df2[tdf1df2[,xveh[i]] <= -500, uveh[i]]
    t     <- tveh[tdf1df2[,xveh[i]] <= -500]
    x.xe[i] <- x[length(t)]
    u.xe[i] <- u[length(t)]
    t.xe[i] <- t[length(t)]
  }
  for(i in 1:10) h[i] <- hsafe(u.x0[i], 14)
  abline(h = c(0,-500), col = gray(0.5))
  abline(v = 0, col = gray(0.5))
  axis(side = 4, at = 0, label = expression(x[0]))
  axis(side = 4, at = -500, label = expression(x[e]))
  df <- data.frame(t.xe, u.xe, t.x0, u.x0, h, x.xe)
  for(i in 1:10) points(t.x0[i], -h[i])
  dt.0 <- t.x0 - t.xe
  df <- cbind(df, dt.0)
  for(i in 2:10) {
    lines(c(t.xe[i],t.x0[i]), c(x.xe[i], x.xe[i] + u.xe[i] * df[i,7]), lwd = 3, col = "red")
  }
  title(main = "SBS Driver Inefficiencies", sub = expression("Decision location: "*x[e]))

### u.2 #####################################################################
  browser()
  u    <- as.ts(data.frame(u1 =  tdf1df2[,2], u2 =  tdf1df2[,5],
                           x1 =  tdf1df2[,3], x2 =  tdf1df2[,6]))
  dt <- 0.125
  Zt <- matrix(c(1,0), 1, 2)
  Ht <- matrix(NA)
  Tt <- matrix(c(1,0,1,1),2,2)
  Rt <- matrix(c(1,0),2,1)
  Qt <- matrix(NA)
  a1 <- matrix(c(1,0),2,1)
  P1 <- matrix(0,2,2)
  P1inf <- diag(2)
  par(mfrow = c(1,2), pty = "s")
  model_u <- SSModel(u2 ~ -1 +
                SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, a1 = a1, P1 = P1, P1inf = P1inf),
                H = Ht, data = u
                )
  fit_u   <- fitSSM(model_u, inits = c(0,0), method = "BFGS")
  out_u   <- KFS(fit_u$model)
  print(out_u)
  print(attributes(out_u))
  plot(u[,1], ylim = c(0, max(c(u[,1],u[,2]))), col = "blue", typ = "p")
  points(u[,2], col = "red")
#  lines(out_u$att[,1], col = "red", lwd = 2)
  lines(out_u$muhat, col = "black", lwd = 2)
  title(main = "Lead Vehicle", sub = "Gaussian")
  legend("bottomleft", legend = c("Observations: 1", "Observations: 2", "Smoothed estimates"),
         lty = c(NA,NA,1), col = c("blue", "red", "black"), lwd = c(NA,NA,2), pch = c(1,1,NA), bty = "n"
  )
  P <- ts(as.matrix(out_u$P), frequency = 1, start = 0, end = 320)
  plot(P, col = "black", typ = "l", lwd = 2)
  title("Error Variance")
  browser()



###########################################################################


}


#' \code{Helske5.model} is a wrapper function for estimating a successful, safe merge.
#'
#' @usage uses \code{Helske5.model}, \code{xab} and \code{uab} from the \code{cartools} Package.
#' @param tdf1df2, time, speed and locations, a matrix
#' @param veh, a number
# #' @examples
# #' Helske5.model(tdf1df2, veh)
#' @export
Helske5.model <- function(tdf1df2, veh) {
  vehseq <- seq(2,31,3)
  z      <- tdf1df2[1:100, vehseq[veh]]
  H      <- round(5 * 5280/3600, 2)
  u      <- 78
  Q0     <- 0
  x      <- rnorm(100, u, Q0)
  df     <- data.frame(x, z)
  start  <- 0
  end    <- 100
  data   <- ts(df, start, end)
  # See Helske3.R
  dt     <- 0.125
  Zt     <- matrix(c(1,dt,dt^2),1,3)
  Ht     <- matrix(NA)
  Tt     <- matrix(c(1,0,0,0,1,0,0,0,1),3,3)
  Rt     <- diag(1,3)
  Qt     <- matrix(c(0,0,0,0,0,0,0,0,NA),3,3)
  model2 <- SSModel(z ~  SSMcustom(Z = Zt, T = Tt, R = Rt, Qt), H = Ht, data = data)
  fit2   <- fitSSM(model2, inits = c(0,0), method = "BFGS")
  out2   <- KFS(fit2$model)
  Qhat2  <- round(as.numeric(fit2$model$Q),2)
  Hhat2  <- round(as.numeric(fit2$model$H),2)

  par(mfrow = c(1,1), pty = "s", mar = c(3,0,1,1))

  ts.plot(data[,1], data[,2], ts(rowSums(out2$a)), out2$muhat,
          col = c("gold", gray(0.5), "blue",  "black"),
          #      ylim = c(u - 3.5*max(c(Q0,H)), u + 3.5*max(c(Q0,H))),
          ylab = "u(t), feet per second", lty = c(1,3,1,1), lwd = c(6,2,2,2)
  )
  title(main = "Intercept Speed Model")
  legend("bottomright",
         legend = c("Gold Standard", "Observed", "One-step ahead", "Smoothed" ),
         lty = c(1,3,1,1),
         lwd = c(6,2,2,2),
         col = c("gold", gray(0.5), "blue","black"),
         bty = "n")

  legend("topleft",c(
    expression(""),
    bquote(u == .(u)),
    bquote(sigma == .(H))),
    bty = "n"
  )
  return(list(model2, fit2, out2, df, out2$v))
}

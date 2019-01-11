#' \code{Helske3.model} is a wrapper function for estimating a successful, safe merge.
#'
#' @usage uses \code{Helske3.model}, \code{xab} and \code{uab} from the \code{cartools} Package.
# #' @param tdf1df2, time, speed and locations, a matrix
# #' @examples
# #' Helske3.model()
#' @export
Helske3.model <- function() {
  Q0     <- 0
  H      <- round(5 * 5280/3600, 2) # Observation noise
  u      <- 78
  x      <- rnorm(100, u, Q0)
  z      <- rnorm(100, x, H)
  df     <- data.frame(x,z)
  start  <- 0
  end    <- 100
  data   <- ts(df, start, end)
  model  <- SSModel(z ~  SSMcustom(Z = matrix(c(1,0.125),1,2),
                                  T = diag(1,2),
                                  R = diag(1,2),
                                  Q = matrix(NA,2,2),
                                  P1 = matrix(0,2,2),
                                  ),
                    H = NA,
                    distribution = "gaussian",
                    data = data)
  updatefn <- function(pars, model) {
    model["H"]  <- pars[1]
    model["Q"]  <- c(pars[2], pars[3])
    model
  }
  fit   <- fitSSM(model,  updatefn = updatefn, inits = as.numeric(c(1,1,1)), method = "BFGS")
  out   <- KFS(fit$model)
  par(mfrow = c(1,1), pty = "s", mar = c(3,0,1,1))
  ts.plot(data[,1], data[,2], ts(rowSums(out$a)), out$muhat,
          col = c("gold", gray(0.5), "blue",  "black"),
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
  if(FALSE) {
    ts.plot(out$P[1,,], col = "black", ylab = "P", lwd = 2)
    title(main = "Covariance")
    legend("topright",c(
      expression(""),
      bquote(H == .(H)),
      bquote(Q == .(Q0)),
      bquote(hat(H) == .(Hhat2)),
      bquote(hat(Q) == .(Qhat2))),
      bty = "n"
    )
  }


  ### Notes
  # 1. u = 50 = speed
  # 2. Q = 0 Here, the goal is to reach u = 50
  # 3. The P array is dificult to interpret.
  # Estimate    Std. Error
  # (Intercept)   7.882e+01   1.897e-08
  # custom1       0.000e+00   0.000e+00
  # custom2      -4.221e-14   1.518e-07

  return(list(model, fit, out, df))
}

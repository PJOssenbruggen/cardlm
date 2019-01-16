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
  # See Helske4.R
  model  <- SSModel(z ~  SSMcustom(Z = matrix(c(1,0.125),1,2),
                                   T = diag(1,2),
                                   R = diag(1,2),
                                   Q = matrix(c(NA,0,0,NA),2,2),
                                   P1 = matrix(c(200,0,0,200),2,2),
            ),
            H = NA,
            distribution = "gaussian",
            data = data,
            tol  = .Machine$double.eps^0.5)
  updatefn <- function(pars, model) {
    model["H"]  <- pars[1]
    model["Q"]  <- c(pars[2], pars[3])
    model
  }
  check_model <- function(model) (model["H"]  > 0 &  model["Q"] > 0)
  fit   <- fitSSM(model,  updatefn = updatefn,
                  check_fn = check_model,
                  inits = as.numeric(c(1,1,1)), method = "BFGS")
  out   <- KFS(fit$model)
  print(out$P[,,end])


  par(mfrow = c(1,1), pty = "s", mar = c(3,0,1,1))

  ts.plot(data[,1], data[,2], ts(rowSums(out$a)), out$muhat,
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

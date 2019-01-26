#' \code{Helske2.model} is a wrapper for estimating speed and location with Kalman filtering \code{SSMcustom} See \code{Helske9.R}.
#'
#' @param tdf1df10 speed and location data, data.frame
# #' @examples
# #' Helske10.model(tdf1df2)
#' @usage Helske10.model(tdf1df2)
#' @export
Helske10.model <- function(tdf1df2) {
  df     <- tdf1df2[,c(2,3,5,6)]
  start  <- 0
  end    <- 40
  u      <- round(53*5280/3600,0)
  usd    <- round(5280/3600*5,1)
  data   <- ts(df, start, end, frequency = 8)
  par(mfrow = c(1,2), pty = "s")
  acf(diff(df[,1]),  main = expression(nabla*dot(x)[t]))

  y  <- ts(diff(df[,1]), start = c(0,1), end = 40, frequency = 8)
  ts.plot(y,
          col = "black",
          ylim = c(-10,10)
  )
  title(main = expression("Arima, "*nabla*dot(x)[t]))
  browser()
  model  <- SSModel(data[,1] ~  -1 +
                      SSMcustom(Z  = matrix(c(1,0.125),1,2),
                                T  = matrix(c(1,0,1,1),2,2),
                                Q  = matrix(NA,1,1),
                                P1 = matrix(c(1,0,0,1),2,2)
                      ),
                    distribution = "gaussian",
                    #           u = data[,3:4],
                    tol  = .Machine$double.eps^0.5,
                    H = NA
  )
  ###########################################################################################################
  check_model  <- function(model) (  model$Q[1,1,1] > 0 &
                                       model$H[1,1,1] > 0
  )
  fit     <- fitSSM(model, inits = c(0,0), checkfn = check_model, method = "BFGS")
  out     <- KFS(fit$model)

  Out <- ts(data.frame(stand = rep(NA,321),
                       obs = data[,1],
                       smooth = signal(out)$signal[[1]],
                       prd    = signal(out,filtered = TRUE)$signal),
            start, end, frequency = 8)
  par(mfrow = c(1,1), pty = "s")
  ts.plot(window(Out, start = c(0,1), end = c(40,0)),
          col = c("gold", gray(0.5), "black", "blue"),
          ylim = c(0,150),
          ylab = "u(t), feet per second", lty = c(1,3,1,1), lwd = c(6,2,2,2)
  )
  title(main = expression(dot(x)[t]))
  legend("bottomright",
         legend = c("target", "Observed", "One-step ahead predictions", "Smoothed estimates" ),
         lty = c(1,3,1,1),
         lwd = c(6,2,2,2),
         col = c("gold", gray(0.5), "blue","black"),
         bty = "n")

  legend("topleft",c(
    expression(""),
    bquote(bar(u) == .(u)),
    bquote(sigma[w] == .(usd))),
    bty = "n"
  )



  ### Notes
  # 1. u = 50 = speed
  # 2. Q = 0 Here, the goal is to reach u = 50
  # 3. tracks well for all usd
  # 4. Covariance quickly reach steady state
  # 5. An acceleration term did not work.
  # 6. No warning message of ldl failure.
  # 7. transform = "augment" halts the analysis.
  # 8. H, a 2x2 matrix, has only one unknown parameters.
  browser()

}

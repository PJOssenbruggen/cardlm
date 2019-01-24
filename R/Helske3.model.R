#' \code{Helske3.model} is a wrapper for estimating speed and location with Kalman filtering \code{SSMcustom}.
#'
#' @usage uses \code{Helske3.model} is used to filter noisy speed data.
#' @param df speed and location data, data.frame
# #' @examples
# #' Helske3.model(df)
#' @export
Helske3.model <- function(df) {
  df     <- df[,1:4]
  start  <- 0
  end    <- 40
  u      <- round(53*5280/3600,0)
  usd    <- round(5280/3600*usd,1)
  data   <- ts(df, start, end, frequency = 8)
  model  <- SSModel(data[,3] ~  -1 +
                      SSMtrend(1, Q = NA) +
                      SSMcustom(Z  = matrix(1,1,1),
                      T  = matrix(1,1,1),
                      Q  = matrix(NA),
                      P1 = diag(1,1,1)
                      ),
                      H = NA,
                      distribution = "gaussian",
                      data = data,
                      tol  = .Machine$double.eps^0.5
                    )

  browser()
  check_model <- function(model) (model["H"]  > 0 &
                                    model["Q", etas = "level"] > 0
                                  )
  fit         <- fitSSM(model,
                        check_fn = check_model,
                        inits    = rep(u/10,3),
                        method   = "BFGS")
  out         <- KFS(fit$model)
  print(out$P[,,end])

  par(mfrow = c(1,2),  pty = "s")
  print("Q and R")
  print(fit$optim.out[[1]])
  print("fit$optim.out")
  print(out)
  print(head(data))
  Out <- ts(data.frame(stand = data[,1], obs = data[,3], smooth = out$muhat, prd = out$att[,1] + 0.125 * out$att[,2]),
            start, end, frequency = 8)

  ts.plot(window(Out, start = c(0,1), end = c(40,0)),
          col = c("gold", gray(0.5), "black", "blue"),
          ylim = c(50,110),
          ylab = "u(t), feet per second", lty = c(1,3,1,1), lwd = c(6,2,2,2)
  )
  title(main = expression(dot(x)[t]))
  legend("bottomright",
         legend = c("Gold Standard", "Observed", "One-step ahead predictions", "Smoothed estimates" ),
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

  p1 <- as.numeric(out$P[1,,][1,])
  p2 <- as.numeric(out$P[1,,][2,])
  tseq <- seq(start, end, 0.125)
  plot(tseq, p1[-1], typ = "l", lwd = 2, ylim = c(0,max(c(p1,p2))), xlab = "Time", ylab = "P")
  lines(tseq,p2[-1], lwd = 2)

  title(main = "Covariance")

  ### Notes
  # 1. u = 50 = speed
  # 2. Q = 0 Here, the goal is to reach u = 50
  # 3. tracks well for all usd
  # 4. Covariance quickly reach steady state
  # 5. Use default updatefn.

  browser()
}

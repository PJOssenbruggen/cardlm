#' \code{Helske3.model} is a wrapper function for estimating a successful, safe merge. Mean speed = 53 mph = 78 fps.
#'
#' @usage uses \code{Helske3.model} is used to filter noisy speed data.
#' @param usd standard deviation of speed in mph, a number.
# #' @examples
# #' Helske3.model()
#' @export
Helske3.model <- function(usd) {
  set.seed(123)
  Q0     <- 0
  H      <- round(usd * 5280/3600,1) # Observation noise
  u      <- 78
  x      <- rnorm(321, u, Q0)        # Gold Standard
  z      <- rnorm(321, x, H)
  df     <- data.frame(x, z)
  start  <- 0
  end    <- 40
  data   <- ts(df, start, end, frequency = 8)

  model  <- SSModel(z ~  -1 +
                      SSMcustom(Z  = matrix(c(1,0.125),1,2),
                      T  = matrix(c(1,0,0,1),2,2),
                      R  = matrix(c(1,0),2,1),
                      Q  = matrix(NA),
                      P1 = diag(1,2,2)
                      ),
                      H = NA,
                      distribution = "gaussian",
                      data = data,
                      tol  = .Machine$double.eps^0.5
                    )

  print("y_t+1 = Z*a_t + H")
  print("Z")
  print(model$Z)
  print("H")
  print(model$H)
  print("a_t+1 = T*a_t + R*Q")
  print("T")
  print(model$T)
  print("R")
  print(model$R)
  print("Q")
  print(model$Q)

  check_model <- function(model) (model["H"]  > 0 &  model["Q"] > 0)
  fit         <- fitSSM(model,
               #         updatefn = updatefn,
                        check_fn = check_model,
                        inits    = rep(u/10,2),
                        method   = "BFGS")
  out         <- KFS(fit$model)
  print(out$P[,,end])

  layout(mat = matrix(c(1,1,2,0),2,2), widths = c(3,1), height = c(3,1))
  par(mar = c(1,3,1,3), pty = "s")
  print("Q and R")
  print(fit$optim.out[[1]])
  print("fit$optim.out")
  print(out)

  Out <- ts(data.frame(stand = data[,1], obs = data[,2], smooth = out$muhat, prd = out$att[,1] + 0.125 * out$att[,2]),
            start, end, frequency = 8)

  ts.plot(window(Out, start = c(0,1), end = c(10,0)),
          col = c("gold", gray(0.5), "black", "blue"),
          ylim = c(u - 3.5*max(c(Q0,H)), u + 3.5*max(c(Q0,H))),
          ylab = "u(t), feet per second", lty = c(1,3,1,1), lwd = c(6,2,2,2)
  )
  title(main = "SSMcustom Speed Model")
  legend("bottomright",
         legend = c("Gold Standard", "Observed", "One-step ahead predictions", "Smoothed estimates" ),
         lty = c(1,3,1,1),
         lwd = c(6,2,2,2),
         col = c("gold", gray(0.5), "blue","black"),
         bty = "n")

  legend("topleft",c(
    expression(""),
    bquote(bar(u) == .(u)),
    bquote(sigma[w] == .(H))),
    bty = "n"
  )
  P  <- ts(out$P[1,,][1,], start, end, frequency = 8)
  ts.plot(window(P, start = c(0,1), end = c(3,0)), col = "black", ylab = "P", lwd = 2)
  title(main = "Covariance")

  ### Notes
  # 1. u = 50 = speed
  # 2. Q = 0 Here, the goal is to reach u = 50
  # 3. tracks well for all usd
  # 4. Covariance quickly reach steady state
  # 5. Use default updatefn.


  return(list(model, fit, out, df))
}

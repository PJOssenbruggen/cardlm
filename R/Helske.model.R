#' \code{Helske.model} is a wrapper function for estimating a successful, safe merge. Mean speed = 53 mph = 78 fps.
#'
#' @param usd standard deviation of speed in mph, a number.
# #' @examples
# #' Helske.model(usd)
#' @usage Helske.model(usd)
#' @export
Helske.model <- function(usd) {
  set.seed(123)
  Q0     <- 0
  H      <- round(usd * 5280/3600,1) # Observation noise
  u      <- 78
  x      <- rnorm(321, u, Q0)        # Gold Standard
  z      <- rnorm(321, x, H)
  df     <- data.frame(x,z)
  start  <- 0
  end    <- 40
  data   <- ts(df, start, end, frequency = 8)

  model  <- SSModel(z ~ SSMtrend(1, Q = NA), H = NA, data = data)
  fit    <- fitSSM(model, inits = c(0,0), method = "BFGS")
  out    <- KFS(fit$model)

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

  layout(mat = matrix(c(1,1,2,0),2,2), widths = c(3,1), height = c(3,1))
  par(mar = c(1,3,1,3), pty = "s")
  ts.plot(data[,1], data[,2], out$att, out$muhat,
          col = c("gold", gray(0.5), "blue",  "black"),
          ylim = c(u - 3.5*max(c(Q0,H)), u + 3.5*max(c(Q0,H))),
          ylab = "u(t), feet per second", lty = c(1,3,1,1), lwd = c(6,2,2,2)
  )
  title(main = "SSMtrend Speed Model")
  legend("bottomright",
         legend = c("Gold Standard", "Observed", "One-step ahead predictions", "Smoothed predictions" ),
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

  ts.plot(ts(out$P[1,,], start, end, frequency = 8), col = "black", ylab = "P", lwd = 2)
  title(main = "Covariance")


  ### Notes
  # 1. u = 53 mph = 78 fps = speed
  # 2. Q = 0 Here, the goal is to reach u = 50
  # 3. tracks well for all usd
  # 4. Covariance quickly reach steady state.

  return(list(model, fit, out, data))
}

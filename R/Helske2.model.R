#' \code{Helske2.model} is a wrapper function for estimating a successful, safe merge. Mean speed = 53 mph = 78 fps.
#'
#' @param usd standard deviation of speed in mph, a number.
# #' @examples
# #' Helske2.model(usd)
#' @usage Helske2.model(usd)
#' @export
Helske2.model <- function(usd) {
  set.seed(123)
  Q0     <- 0
  H      <- round(usd * 5280/3600,1) # Observation noise
  u      <- 78
  tseq   <- seq(0,40,0.125)
  u.     <- rnorm(321, u, Q0)        # Gold Standard
  x.     <- u * tseq - 700
  u.r    <- rnorm(321, u, H)
  x.r    <- rep(NA, length(tseq))
  x.r[1] <- -700
  for(i in 2:length(tseq)) x.r[i] <- x.r[i-1] + u.r[i] * 0.125
  df     <- data.frame(u.gold = u., x.gold = x., z1 = u.r, z2 = x.r)

  start  <- 0
  end    <- 40
  data   <- ts(df, start, end, frequency = 8)

  model  <- SSModel(data[,c(3,4)] ~  -1 +
                      SSMcustom(Z  = matrix(c(1,0.125,0,1),2,2),
                                T  = matrix(c(1,0.125,0,1),2,2),
                                R  = matrix(c(1,0.125/2),2,1),
                                Q  = matrix(NA),
                                P1inf = diag(2)
        ),
        H = matrix(c(NA,0,0,0),2,2),
        distribution = "gaussian",
        data = data,
        tol  = .Machine$double.eps^0.5)
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
###########################################################################################################

  check_model  <- function(model) (model$H[1,1,1] > 0 &  model["Q"] > 0)
  update_model <- function(pars, model) {
    model$H[1,1,1] <- pars[1]
    model["Q"]     <- pars[2]
    model
  }
  fit         <- fitSSM(model,
                        updatefn = update_model,
                        checkfn  = check_model,
                        inits    = rep(10,2),
                        method   = "BFGS")
  out         <- KFS(fit$model, transform = "augment")
  print(out$P[,,end])
  print(out$Pinf)

  layout(mat = matrix(c(1,1,2,3),2,2), widths = c(3,1), height = c(3,1))
  par(mar = c(1,3,1,3), pty = "s")
  Out <- ts(data.frame(df, smooth = out$muhat, prd = out$att), start, end, frequency = 8)
  ts.plot(window(Out, start = c(0,1), end = c(10,0))[,c(1,3,5,7)],
          col = c("gold", gray(0.5), "black", "blue"),
          ylim = c(u - 3.5*max(c(Q0,H)), u + 3.5*max(c(Q0,H))),
          ylab = "u(t), feet per second", lty = c(1,3,1,1), lwd = c(6,2,2,2)
  )
  title(main = "SSMcustom Speed Model Predictions")
  legend("bottomright",
         legend = c("Gold Standard", "Observed", "One-step ahead", "Smoothed" ),
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
  ts.plot(window(Out, start = c(0,1), end = c(30,0))[,c(2,4,6,8)],
          col = c("gold", gray(0.5), "black", "blue"),
          ylim = c(-800,1500),
          ylab = "x(t), feet", lty = c(1,3,1,1), lwd = c(6,2,2,2)
  )
  title(main = "Location Predictions")
  Pout <- data.frame(out$P[2,,][1,], out$P[2,,][2,])
  P  <- ts(Pout, start = c(0,1), end = c(10,1), frequency = 8)
  ts.plot(window(P, start = c(0,1), end = c(3,0)), col = "black", ylab = "P", lwd = 2)
  title(main = "Covariance")

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
  return(list(model, fit, out, df))
}

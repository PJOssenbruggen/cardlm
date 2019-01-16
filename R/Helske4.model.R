#' \code{Helske4.model} is a wrapper function for estimating a successful, safe merge.
#'
#' @param usd standard deviation of speed in mph, a number.
#' @usage uses \code{Helske4.model}, \code{xab} and \code{uab} from the \code{cartools} Package.
#' @param veh, a number
# #' @examples
# #' Helske4.model(usd)
#' @export
Helske4.model <- function(usd) {
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

  model  <- SSModel(data[,c(3,4)] ~ SSMtrend(1, Q = matrix(NA), index = c(1,2)),
                      H    = matrix(c(NA,0,0,NA),2,2),
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

  browser()

  check_model  <- function(model) (model["H"]  > 0 &  model["Q"] > 0)
  update_model <- function(pars, model) {
       model["H"] <- pars[1:2]
       model["Q"] <- pars[3:6]
       model
       }
  fit         <- fitSSM(model,
                        updatefn = update_model,
                        checkfn  = check_model,
                        inits    = rep(0.1,6),
                        method   = "BFGS")
  out         <- KFS(fit$model)
  print(out$P[,,end])

  layout(mat = matrix(c(1,1,2,0),2,2), widths = c(3,1), height = c(3,1))
  par(mar = c(1,3,1,3), pty = "s")
  print("Q and R")
  print(fit$optim.out[[1]])

  browser()

  Out <- ts(data.frame(stand = data[,1], obs = data[,3], smooth = out$muhat[,1], prd = out$a[-1,1]),
            start, end, frequency = 8)
  ts.plot(window(Out, start = c(0,1), end = c(10,0)),
          col = c("gold", gray(0.5), "black", "blue"),
          ylim = c(u - 3.5*max(c(Q0,H)), u + 3.5*max(c(Q0,H))),
          ylab = "u(t), feet per second", lty = c(1,3,1,1), lwd = c(6,2,2,2)
  )
  title(main = "Location/Speed Trend Model")
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
  P  <- ts(out$P[1,,][1,], start, end, frequency = 8)
  ts.plot(window(P, start = c(0,1), end = c(3,0)), col = "black", ylab = "P", lwd = 2)
  title(main = "Covariance")

  ### Notes
  # 1. u = 50 = speed
  # 2. Q = 0 Here, the goal is to reach u = 50
  # 3. tracks well for all usd
  # 4. Covariance quickly reach steady state


  return(list(model, fit, out, df))
}

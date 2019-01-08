#' \code{Helske.model} is a wrapper function for estimating a successful, safe merge.
#'
#' @usage uses \code{Helske.model}, \code{xab} and \code{uab} from the \code{cartools} Package.
# #' @param tdf1df2, time, speed and locations, a matrix
# #' @examples
# #' Helske.model()
#' @export
Helske.model <- function() {
  set.seed(1)
  x     <- cumsum(rnorm(100, 0, 0.1))  # No noise
  y     <- rnorm(100, x, 2)          # Noise added
  df    <- data.frame(x,y)

  start <- 0
  end   <- 100
  data  <- ts(df, start, end)
  data

  model <- SSModel(y ~ SSMtrend(1, Q = 0.01), H = 0.01)
  out   <- KFS(model)
  par(mfrow = c(1,1), pty = "s")
  ts.plot(ts(x), out$a, out$att, out$alpha, col = 1:4,  ylab = "y")
  title(main = "out")
  browser()
  print(data.frame(Z = as.numeric(model$Z), H = as.numeric(model$H), T = as.numeric(model$T),
            Q = as.numeric(model$Q)))
  model1 <- SSModel(y ~ SSMtrend(1, Q = NA), H = NA)
  fit1   <- fitSSM(model1, inits = c(0,0), method = "BFGS")
  out1   <- KFS(fit1$model)
  print(data.frame(Z = as.numeric(fit1$model$Z), H = as.numeric(fit1$model$H), T = as.numeric(fit1$model$T),
                   Q = as.numeric(fit1$model$Q)))


  simmodel <- SSModel(y ~ SSMtrend(1, Q = as.numeric(fit1$model$Q)), H = as.numeric(fit1$model$H))
  sim      <- simulateSSM(simmodel, type = "states")

  ts.plot(ts(x), out1$a, out1$att, out1$alpha, ts(sim[,,1]),
          col = c(gray(0.5), "red", "blue", "wheat", "black"),
          ylab = "y", lty = c(1,1,1,1,2), lwd = c(1,1,1,1,4))
  title(main = "out1")
  print(data.frame(ts(x), out1$a[-1], out1$att, out1$alpha,ts(sim[,,1])))
  browser()

  ts.plot(out$P[1,,], col = "black", ylab = "P")
  title(main = "out")

  ts.plot(out1$P[1,,], col = "black", ylab = "P")
  title(main = "out1")
  browser()
#  model2 <- SSModel(y ~ SSMtrend(1, Q = NA),
#          SSMcustom(Z = 1, T = 0, Q = NA, P1 = NA),
#          distribution = "poisson", u = x)
#  update_poisson <- function(pars, model2) {
#    model2["Q", etas = "level"] <- exp(pars[1])
#    model2["Q", etas = "custom"] <- exp(pars[1])
#    model2
#  }
#  fit2   <- fitSSM(model2, inits = c(3,3), updatafn = update_poisson, method = "BFGS")
#  fit2$model["Q", etas = "level"]
#  fit2$model["Q", etas = "custom"]
#  out2   <- KFS(fit2$model)
#  ts.plot(ts(x), out2$a, out2$att, out2$alpha, col = 1:4)
#  title(main = "out2")
  return(list(model1, fit1, out1))
}

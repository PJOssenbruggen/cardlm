#' \code{Helske.model} is a wrapper function for estimating a successful, safe merge.
#'
#' @usage uses \code{Helske.model}, \code{xab} and \code{uab} from the \code{cartools} Package.
# #' @param tdf1df2, time, speed and locations, a matrix
# #' @examples
# #' Helske.model()
#' @export
Helske.model <- function() {
  Q0     <- 0
  H      <- round(5 * 5280/3600,1) # Observation noise
  u      <- 78
  x      <- rnorm(100, u, Q0)
  z      <- rnorm(100, x, H)
  df     <- data.frame(x,z)
  start  <- 0
  end    <- 100
  data   <- ts(df, start, end)
  model1 <- SSModel(z ~ SSMtrend(1, Q = NA), H = NA, data = data)
  fit1   <- fitSSM(model1, inits = c(0,0), method = "BFGS")
  out1   <- KFS(fit1$model)
  Qhat1  <- round(as.numeric(fit1$model$Q),2)
  Hhat1  <- round(as.numeric(fit1$model$H),2)
  df     <- data.frame(Z = as.numeric(fit1$model$Z), H, Hhat1, Q = Q0, Qhat1)

  #par(mfrow = c(1,1), pty = "s")
  # plot 1
  layout(mat = matrix(c(1,1,0,2),2,2), widths = c(3,1), height = c(3,1))
  par(mar = c(3,0,1,1))
  ts.plot(data[,1], data[,2], out1$a, out1$muhat,
          col = c("gold", gray(0.5), "blue",  "black"),
          ylim = c(u - 3.5*max(c(Q0,H)), u + 3.5*max(c(Q0,H))),
          ylab = "u(t), feet per second", lty = c(1,3,1,1), lwd = c(6,2,2,2)
  )
  title(main = "Predictions")
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

  ts.plot(out1$P[1,,], col = "black", ylab = "P", lwd = 2)
  title(main = "Covariance")
  legend("topright",c(
    expression(""),
    bquote(H == .(H)),
    bquote(Q == .(Q0)),
    bquote(hat(H) == .(Hhat1)),
    bquote(hat(Q) == .(Qhat1))),
    bty = "n"
  )


  ### Notes
  # 1. u = 50 = speed
  # 2. Q = 0 Here, the goal is to reach u = 50

  if(FALSE) {
    ### Q and H assigned
    model <- SSModel(z ~ SSMtrend(1, Q = 0.01), H = 0.01)
    out   <- KFS(model)
    par(mfrow = c(1,1), pty = "s")
    ts.plot(ts(x), out$a, out$att, out$alpha, col = 1:4,  ylab = "z")
    title(main = "out")
    browser()
    print(data.frame(Z = as.numeric(model$Z), H = as.numeric(model$H), T = as.numeric(model$T),
                     Q = as.numeric(model$Q)))
    ### Cannot debug the following "poisson."
    ts.plot(out1$P[1,,], col = "black", ylab = "P")
    title(main = "out1")

    model2 <- SSModel(z ~ SSMtrend(1, Q = NA),
                      SSMcustom(Z = 1, T = 0, Q = NA, P1 = NA),
                      distribution = "poisson", u = x)
    update_poisson <- function(pars, model2) {
      model2["Q", etas = "level"] <- exp(pars[1])
      model2["Q", etas = "custom"] <- exp(pars[1])
      model2
    }
    fit2   <- fitSSM(model2, inits = c(3,3), updatafn = update_poisson, method = "BFGS")
    fit2$model["Q", etas = "level"]
    fit2$model["Q", etas = "custom"]
    out2   <- KFS(fit2$model)
    ts.plot(ts(x), out2$a, out2$att, out2$alpha, col = 1:4)
    title(main = "out2")
    # 3. simulate
    simmodel <- SSModel(z ~ SSMtrend(1, Q = as.numeric(fit1$model$Q)), H = as.numeric(fit1$model$H))
    sim      <- simulateSSM(simmodel, type = "states")
  }

  return(list(model1, fit1, out1, df, out1$v))
}

#' \code{Helske2.model} is a wrapper for estimating speed and location with Kalman filtering \code{SSMcustom}.
#'
#' @param df speed and location data, data.frame
# #' @examples
# #' Helske2.model(df)
#' @usage Helske2.model(df)
#' @export
Helske2.model <- function(df) {

  df     <- df[,1:4]
  start  <- 0
  end    <- 40
  u      <- 78
  data   <- ts(df, start, end, frequency = 8)
  model  <- SSModel(data[,c(3,4)] ~  -1 +
                      SSMtrend(2, Q = list(matrix(NA,2,2), matrix(0,2,2))) +
                      SSMcustom(Z  = matrix(c(1,0,1,0.125),2,2),
                                T  = matrix(c(1,0,0,1),2,2),
                                Q  = matrix(NA,2,2),
                                P1 = matrix(NA,2,2)
                      ),
                    distribution = "gaussian",
                    u = data[,3:4],
                    tol  = .Machine$double.eps^0.5
        )
###########################################################################################################
  check_model  <- function(model) (  model$Q[1,1,1] > 0 &
                                     model$Q[2,2,1] > 0 &
                                     model$Q[3,3,1] > 0
                                   )
  update_model <- function(pars, model) {
    browser()
    Q               <- diag(pars[1:2])
    Q[upper.tri(Q)] <- 0                                   # ?????????????????
    model["Q", etas = "level"] <- crossprod(Q)
    Q               <- pars[3]
    model["Q", etas = "custom"] <- model["P1", states = "custom"] <- Q
    print(model["Q"])
    model
  }
  browser()
  init     <- chol(cov(data[,3:4]))
  fitinit  <- fitSSM(model,
                     updatefn = update_model,
                     inits    = rep(c(diag(init), init[upper.tri(init)]),1),
                     method   = "BFGS")
  print(-fitinit$optim.out$val)
  fit <- fitSSM(model,
                updatefn = update_model,
                inits = fitinit$optim.out$par,
                method = "BFGS", nsim = 250)
  print(-fitinit$optim.out$val)
  varcor <- fit$model["Q", etas = "level"]
  varcor[upper.tri(varcor)] <- cov2cor(varcor)[upper.tri(varcor)]
  print(varcor, digits = 2)

  varcor <- fit$model["Q", etas = "custom"]
  varcor[upper.tri(varcor)] <- cov2cor(varcor)[upper.tri(varcor)]
  print(varcor, digits = 2)
  out <- KFS(fit$model, nsim = 1000)
  print(out)

  browser()
  out         <- KFS(fit$model, transform = "augment")
  print(out$P[,,end])
  print(out$Pinf)

  layout(mat = matrix(c(1,1,2,3),2,2), widths = c(3,1), height = c(3,1))
  par(mar = c(1,3,1,3), pty = "s")
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
  p3 <- as.numeric(out$P[1,,][3,])
  plot(tseq, p1[-1], typ = "l", lwd = 2, ylim = c(0,max(c(p1,p2,p3))), xlab = "Time", ylab = "P")
  lines(tseq,p2[-1], lwd = 2)
  lines(tseq,p3[-1], lwd = 2)
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

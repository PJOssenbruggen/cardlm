#' \code{Helske.model} is a wrapper for estimating speed and location with Kalman filtering \code{SSMtrend}.
#'
#' @param df speed and location data, data.frame
# #' @examples
# #' Helske.model(df)
#' @usage Helske.model(df)
#' @export
Helske.model <- function(df) {
  df     <- df[,1:4]
  start  <- 0
  end    <- 40
  u      <- round(53*5280/3600,0)
  usd    <- round(5280/3600*usd,1)
  data   <- ts(df, start, end, frequency = 8)
  model  <- SSModel(data[,3] ~ SSMtrend(1, Q = NA), H = NA, data = data)
  fit    <- fitSSM(model, inits = c(0,0), method = "BFGS")
  out    <- KFS(fit$model)
  print(out)
  Q      <- fit$optim.out[[1]][1]
  H      <- fit$optim.out[[1]][2]
  df.est <- data.frame(Q,H)
  model.cal <- SSModel(data[,3] ~ SSMtrend(1, Q = fit$optim.out[[1]][1]), H = fit$optim.out[[1]][2], data = data)
  predict(model.cal, interval = "prediction", se.fit = TRUE, level = 0.95, n.ahead = 2)
 ######################################################################################
  model2  <- SSModel(data[,4] ~ SSMtrend(1, Q = NA), H = NA, data = data)
  fit2    <- fitSSM(model2, inits = c(0,0), method = "BFGS")
  out2    <- KFS(fit2$model)
  print(out2)
  Q       <- fit2$optim.out[[1]][1]
  H       <- fit2$optim.out[[1]][2]
  df.est  <- rbind(df.est, data.frame(Q,H))
  print(df.est)

########################################################################################


#########################################################################################
  browser()
  Out <- ts(data.frame(gold.v   = df[,1],
                       gold.z   = df[,2],
                       obs.v    = df[,3],
                       obs.z    = df[,4],
                       smooth.v = out$att,
                       smooth.z = out2$att,
                       prd.v    = out$muhat,
                       prd.x    = out2$muhat),
            start, end, frequency = 8)
  print(head(Out))

  layout(mat = matrix(c(1,1,2,3),2,2), widths = c(3,3), height = c(3,3))

  ts.plot(window(Out, start = c(0,1), end = c(40,1))[,c(1,3,5,7)],
          col = c("gold", gray(0.5), "black", "blue"),
          ylab = "u(t), feet per second", lty = c(1,3,1,1), lwd = c(6,2,2,2)
  )

  title(main = expression(dot(x)[t]))
  legend("bottomright",
         legend = c("Gold Standard", "Observed", "One-step ahead predictions", "Smoothed predictions" ),
         lty = c(1,3,1,1),
         lwd = c(6,2,2,2),
         col = c("gold", gray(0.5), "blue","black"),
         bty = "n")

  legend("topleft",c(
    expression(""),
    bquote(bar(u) == .(u)),
    bquote(sigma[epsilon] == .(usd))),
    bty = "n"
  )


  ts.plot(window(Out, start = c(0,1), end = c(2,1))[,c(2,4,6,8)],
          col = c("gold", gray(0.5), "black", "blue"),
          ylab = "x(t), feet", lty = c(1,3,1,1), lwd = c(6,2,2,2)
  )
  title(main = expression(x[t]))
  legend("topleft",
         legend = c("Gold Standard", "Observed", "One-step ahead predictions", "Smoothed predictions" ),
         lty = c(1,3,1,1),
         lwd = c(6,2,2,2),
         col = c("gold", gray(0.5), "blue","black"),
         bty = "n")

  Out2 <- ts(data.frame(gold.v   = out$P[1,,], gold.z   = out2$P[1,,]),
                     start, end, frequency = 8)

  ts.plot(window(Out2, start = c(0,1), end = c(40,1))[,c(1,2)],
          col = c("black", "blue"),
          ylab = "P", lty = c(1,1), lwd = c(2,2)
  )
  title(main = "Covariance")
  legend("right",
         legend = c("Speed", "Location" ),
         lty = c(1,1),
         lwd = c(2,2),
         col = c("black","blue"),
         bty = "n")

  ### Notes
  # 1. u = 53 mph = 78 fps = speed
  # 2. Q = 0 Here, the goal is to reach u = 50
  # 3. tracks well for all usd
  # 4. Covariance quickly reach steady state.

  browser()

}

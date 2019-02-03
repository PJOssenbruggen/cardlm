#' \code{Helske.model} is a wrapper for estimating speed and location with Kalman filtering \code{SSMtrend}.
#'
#' @param df speed and location data derived from \code{RingRoaddata}, data.frame
# #' @examples
# #' Helske.model(df)
#' @usage
#'
#' @export
Helske.model <- function(df) {
  df     <- df[,1:4]
  start  <- 0
  end    <- 40
  u      <- round(53*5280/3600,0)
  usd    <- round(5280/3600*usd,1)
  data   <- ts(df, start, end, frequency = 8)
  par(mfrow = c(1,2), pty = "s")
  model  <- SSModel(window(data[,3], start = c(0,1), end = c(39,0)) ~
                      SSMtrend(1, Q = NA), H = NA, data = data)
  update_model <- function(pars, model) {
    model["H"] <- pars[1]
    model["Q"] <- pars[2]
    model
  }
  check_model  <- function(model) (model["H"] > 0 &
                                     model$Q[1,1,1] > 0)
  fit         <- fitSSM(model,
                        inits    = rep(0.1,3),
                        checkfn  = check_model,
                        updatefn = update_model,
                        method   = "BFGS")
  out         <- KFS(fit$model)
  Q           <- fit$optim.out[[1]][1]
  H           <- fit$optim.out[[1]][2]
  df.est      <- data.frame(Q,H)
  print(df.est)
  pred  <- predict(fit$model, interval = "prediction",
                    se.fit = TRUE, level = 0.95, n.ahead = 8)
  Out   <- cbind(gold.v = data[,1],
                       obs.v    = data[,3],
                       smooth.v = signal(out,filtered = TRUE)[[1]],
                       prd.v    = signal(out)[[1]],
                       pred     = pred[,1:3]
                      )
  ts.plot(window(Out, start = c(0,1), end = c(40,1)),
          col = c("gold", gray(0.5), "black", "blue", "red","red","red"),
          ylab = "u(t), feet per second",
          lty = c(1,3,1,1,1,2,2),
          lwd = c(6,2,2,2,4,2,2)
  )
  title(main = expression(dot(x)[t]))
  legend("bottom",
         legend = c("Target", "Observed", "Filtered", "Signal","Prediction","95% Prediction interval" ),
         lty = c(1,3,1,1,1,3),
         lwd = c(6,2,2,2,4,2,2),
         col = c("gold", gray(0.5), "blue","black","red","red","red"),
         bty = "n")

  legend("topleft",c(
    expression(""),
    bquote(bar(u) == .(u)),
    bquote(sigma[U] == .(usd))),
    bty = "n"
  )
  browser()
 ######################################################################################

  model  <- SSModel(window(data[,4], start = c(0,1), end = c(39,1)) ~
                           SSMtrend(1, Q = NA), H = NA, data = data)
  update_model <- function(pars, model) {
  model["H"] <- pars[1]
  model["Q"] <- pars[2]
  model
  }
  check_model  <- function(model) (model["H"] > 0 &
                                   model$Q[1,1,1] > 0)
  fit         <- fitSSM(model,
                      inits    = rep(0.1,3),
                      checkfn  = check_model,
                      updatefn = update_model,
                      method   = "BFGS")
  out         <- KFS(fit$model)
  Q           <- fit$optim.out[[1]][1]
  H           <- fit$optim.out[[1]][2]
  df.est      <- data.frame(Q,H)
  print(df.est)
  pred  <- predict(fit$model, interval = "prediction",
                 se.fit = TRUE, level = 0.95, n.ahead = 8)

  Out   <- cbind(gold.v = data[,2],
               obs.v    = data[,4],
               smooth.v = signal(out,filtered = TRUE)[[1]],
               prd.v    = signal(out)[[1]],
               pred     = pred[,1:3]
  )

  ts.plot(window(Out, start = c(39,1), end = c(40,1)),
          ylim = c(2900,3200),
          col = c("gold", gray(0.5), "black", "blue", "red","red","red"),
          ylab = "u(t), feet per second",
          lty = c(1,3,1,1,1,2,2),
          lwd = c(6,2,2,2,4,2,2)
  )
  title(main = expression(x[t]))
browser()
}

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
  usd    <- round(5280/3600*5,1)
  target <- rep(u,dim(df)[1])
  df     <- cbind(df, target)
  data   <- ts(df, start, end, frequency = 8)
  par(mfrow = c(1,2), pty = "s")
  ### Data Set 1 ###################################################################################
  model  <- SSModel(window(data[,1], start = c(0,1), end = c(25,0)) ~ -1 +
                      SSMcustom(Z = matrix(c(1,0), 1, 2),
                                R = matrix(c(1,0,0,1), 2, 2),
                                T = matrix(c(1,0,0,0), 2, 2),
                                Q = matrix(c(NA,0,0,0), 2, 2),
                                a1 = matrix(c(0,0),2,1),
                                P1inf = diag(2),
                                P1 = diag(log(var(data[,1])), 2)),

                    H    = matrix(NA,1,1),
                    data = data
  )

  update_model <- function(pars, model) {
    model$Q[1,1,1]    <- exp(pars[1])
    model$H[1,1,1]    <- exp(pars[2])
    model
  }
  check_model  <- function(model) (model$H[1,1,1] > 0 & model$Q[1,1,1] > 0)
  fit         <- fitSSM(model,
                        inits    = c(log(var(data[,1])),log(var(data[,1]))),
                        checkfn  = check_model,
                        updatefn = update_model,
                        method   = "BFGS")
  out         <- KFS(fit$model)
  Q           <- fit$optim.out[[1]][1]
  H           <- fit$optim.out[[1]][2]
  Q           <- round(Q,1)
  H           <- round(H,1)
  df.est      <- data.frame(Q,H)
  print(df.est)

  pred  <- predict(fit$model, interval = "prediction",
                   se.fit = TRUE, level = 0.95, n.ahead = 160)
  Out   <- cbind(gold.v = data[,5],
                       obs.v    = data[,1],
                       smooth.v = signal(out,filtered = TRUE)[[1]],
                       prd.v    = signal(out)[[1]],
                       pred     = pred[,1:3]
                      )
  ts.plot(window(Out, start = c(0,1), end = c(40,1)),
          ylim = c(0,160),
          col = c("gold", gray(0.5), "black", "blue", "red","red","red"),
          ylab = "u(t), feet per second",
          lty = c(1,3,1,1,1,2,2),
          lwd = c(6,2,2,2,4,2,2)
  )

  title(main = expression(dot(x)[t]), sub = expression("speed/acceleration"))
  legend("bottom",
         legend = c("Target", "Observed", "Filtered", "Signal","Prediction","95% Prediction interval" ),
         lty = c(1,3,1,1,1,3),
         lwd = c(6,2,2,2,4,2,2),
         col = c("gold", gray(0.5), "blue","black","red","red","red"),
         bty = "n")
  legend("topleft",c(
    expression(""),
    bquote(bar(u) == .(u)),
    bquote(sigma[U] == .(usd)),
    bquote(hat(Q) == .(Q)),
    bquote(hat(H) == .(H))
  ),
  bty = "n"
  )
  browser()
 ######################################################################################
  model  <- SSModel(window(data[,1], start = c(0,1), end = c(35,0)) ~ -1 +
                      SSMcustom(Z = matrix(1,1,1),
                                R = matrix(1,1,1),
                                T = matrix(1,1,1),
                                Q = matrix(NA),
                                a1 = matrix(1,1,1),
                                P1inf = diag(1),
                                P1 = diag(log(var(data[,1])), 1)),
                    H    = matrix(NA,1,1),
                    data = data
  )

  update_model <- function(pars, model) {
    model["Q"] <- exp(pars[1])
    model["H"] <- exp(pars[2])
    model
  }
  check_model  <- function(model) (model["H"] > 0 & model["Q"] > 0)
  fit         <- fitSSM(model,
                        inits    = c(log(var(data[,1])),log(var(data[,1]))),
                        checkfn  = check_model,
                        updatefn = update_model,
                        method   = "BFGS")
  out         <- KFS(fit$model)
  Q           <- fit$optim.out[[1]][1]
  H           <- fit$optim.out[[1]][2]
  Q           <- round(Q,1)
  df.est      <- data.frame(Q,H)
  print(df.est)
  pred  <- predict(fit$model, interval = "prediction",
                   se.fit = TRUE, level = 0.95, n.ahead = 80)

  Out   <- cbind(gold.v   = data[,5],
                 obs.v    = data[,1],
                 smooth.v = signal(out, states = "all", filtered = TRUE)[[1]],
                 prd.v    = signal(out, states = "all")[[1]],
                 pred     = pred[,1:3]
  )

  ts.plot(window(Out, start = c(0,1), end = c(45,1)),
          ylim = c(0,150),
          col = c("gold", gray(0.5), "black", "blue", "red","red","red"),
          ylab = "u(t), feet per second",
          lty = c(1,3,1,1,1,2,2),
          lwd = c(6,2,2,2,4,2,2)
  )
  title(main = expression(dot(x)[t]), sub = "SSMcustom")
  legend("topleft",c(
    expression(""),
    bquote(bar(u) == .(u)),
    bquote(sigma[U] == .(usd)),
    bquote(hat(Q) == .(Q)),
    bquote(hat(H) == .(H))
  ),
  bty = "n"
  )
  browser()
}

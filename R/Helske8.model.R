
#' \code{Helske8.model} is a wrapper for learning how to fit \code{SSMarima} models.
#'
#' @param tdf1df2 speed and location data, data.frame
# #' @examples
# #' Helske8.model(tdf1df2)
#' @export
Helske8.model <- function(tdf1df2) {
  df     <- tdf1df2[,c(2,3,5,6)]
  start  <- 0
  end    <- 40
  u      <- round(53*5280/3600,0)
  usd    <- round(5280/3600*5,1)
  data   <- ts(df, start, end, frequency = 8)

  ### SSMcustom Finland_____________________________________________________
  if(FALSE) {
    data("alcohol", package = "KFAS")
    deaths <- window(alcohol[, 2], end = 2007)
    population <- window(alcohol[, 6], end = 2007)
    y     <- deaths/population
    u     <- rep(NA,length(deaths))
    par(mfrow = c(1,2), pty = "s")
    model <- SSModel(deaths/population ~ -1 + SSMcustom(
      Z = matrix(c(1,0),1,2),
      T = matrix(c(1,0,1,1),2,2),
      R = matrix(c(1,0),2,1),
      Q = matrix(NA),
      a1 = matrix(c(1,0),1,2),
      P1 = matrix(0,2,2),
      P1inf = diag(2)
      ),
      H = NA)
    fit   <- fitSSM(model, inits = c(0,0), method = "BFGS")
    out   <- KFS(fit$model)
    Out   <- ts(data.frame(u, obs = y,
                         smooth = signal(out)$signal,
                         prd = signal(out,filtered = TRUE)$signal
                         ),
              start = 1969, end = 2007, frequency = 1)

    ts.plot(window(Out, start = 1969, end = 2007),
            col = c("yellow",gray(0.5), "black", "blue"),
            ylim = c(0,60),
            ylab = "u(t), feet per second", lty = c(1,3,1,1), lwd = c(6,2,2,2)
    )
    title(main = expression("death/population, "*dot(x)[t]))
    legend("bottomright",
           legend = c("Target", "Observed", "One-step ahead predictions", "Smoothed estimates" ),
           lty = c(1,3,1,1),
           lwd = c(6,2,2,2),
           col = c("gold", gray(0.5), "blue","black"),
           bty = "n")
    print(model$T)
    print(model$R)
    print(model$Q)
    browser()
  }
  ### SSMarima Finland______________________________________________________
  if(FALSE) {
    data("alcohol", package = "KFAS")
    deaths <- window(alcohol[, 2], end = 2007)
    population <- window(alcohol[, 6], end = 2007)
    y     <- deaths/population
    u     <- rep(NA,length(deaths))
    drift <- 1:length(deaths)
    model_arima <- SSModel(deaths / population ~ drift +
                                + SSMarima(ma = 0, d = 1, Q = 1))
    update_model <- function(pars, model) {
          tmp <- SSMarima(ma = pars[1], d = 1, Q = pars[2])
          model["R", states = "arima"]  <- tmp$R
          model["Q", states = "arima"]  <- tmp$Q
          model["P1", states = "arima"] <- tmp$P1
          model
      }
    fit_arima <- fitSSM(model_arima, inits = c(0, 1),
                        updatefn = update_model,
                        method = "L-BFGS-B", lower = c(-1, 0), upper = c(1, 100))
    fit_arima$optim.out$par
    out_arima <- KFS(fit_arima$model)
    out_arima$logLik
    Out   <- ts(data.frame(u, obs = y,
                           smooth = signal(out_arima)$signal,
                           prd = signal(out_arima,filtered = TRUE)$signal
                          ),
                          start = 1969, end = 2007, frequency = 1)
    ts.plot(window(Out, start = 1969, end = 2007),
            col = c("yellow",gray(0.5), "black", "blue"),
            ylim = c(0,60),
            ylab = "u(t), feet per second", lty = c(1,3,1,1), lwd = c(6,2,2,2)
    )
    title(main = expression("Arima, "*dot(x)[t]))
    legend("bottomright",
           legend = c("Target", "Observed", "One-step ahead predictions", "Smoothed estimates" ),
           lty = c(1,3,1,1),
           lwd = c(6,2,2,2),
           col = c("gold", gray(0.5), "blue","black"),
           bty = "n")
    browser()
  }
  ### SSMarima y simulated__________________________________________________
  if(TRUE) {
    y      <- arima.sim(n = 321, list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)),
                        innov = rnorm(1000) * sqrt(0.5))
    y      <- ts(as.numeric(y), start, end, frequency = 8)
    u      <- mean(y)
    par(mfrow = c(2,2), pty = "s")
    acf(y)
    drift <- 1:length(y)
    model_arima <- SSModel(y ~ drift +
                             SSMarima(ar = c(0, 0), d = 0, ma = c(0, 0), Q = 1)
                           )
    update_model <- function(pars, model) {
      tmp <- SSMarima(ma = artransform(pars[1:2]),
                      d = 0,
                      ar = artransform(pars[3:4]),
                      Q = exp(pars[5]))
      model["T", states = "arima"] <- tmp$T
      model["R", states = "arima"] <- tmp$R
      model["Q", states = "arima"] <- tmp$Q
      model["P1",states = "arima"] <- tmp$P1
      model
    }
    fit_arima <- fitSSM(model_arima, inits = rep(1,5), updatefn = update_model,
                        method = "L-BFGS-B", lower = c(-1, 0), upper = c(1, 100))
    print(fit_arima$optim.out$par)
    out_arima <- KFS(fit_arima$model)
    Out   <- ts(data.frame(u, obs = y,
                           smooth = signal(out_arima)$signal,
                           prd = signal(out_arima,filtered = TRUE)$signal
              ), start = 1969, end = 2007, frequency = 1)
    ts.plot(window(Out, start = 1969, end = 2007),
            col = c("yellow",gray(0.5), "black", "blue"),
            ylim = c(-20,10),
            ylab = "u(t), feet per second", lty = c(1,3,1,1), lwd = c(6,2,2,2)
    )
    title(main = expression("Arima, "*dot(x)[t]))
    legend("bottomright",
           legend = c("Target", "Observed", "One-step ahead predictions", "Smoothed estimates" ),
           lty = c(1,3,1,1),
           lwd = c(6,2,2,2),
           col = c("gold", gray(0.5), "blue","black"),
           bty = "n")
    browser()
  }
  ### SSMarima y simulated acf agree _______________________________________
  if(TRUE) {
    y      <- arima.sim(n = 321, list(ar = c(0.8897, -0.4858), ma = c(-0.2279, 0.2488)),
                        innov = rnorm(1000) * sqrt(0.5))
    y      <- ts(as.numeric(y), start, end, frequency = 8)
    u      <- mean(y)

    acf(y)
    drift <- 1:length(y)
    model_arima <- SSModel(y ~ drift +
                             SSMarima(ar = c(0, 0), d = 0, ma = c(0, 0), Q = 1)
    )
    likfn <- function(pars, model, estimate = TRUE){
      tmp <- try(SSMarima(artransform(pars[1:2]), d = 0,
                          artransform(pars[3:4]),
                          Q = exp(pars[5])), silent = TRUE)
      if(!inherits(tmp, "try-error")){
        model["T", "arima"]  <- tmp$T
        model["R", "arima"]  <- tmp$R
        model["P1", "arima"] <- tmp$P1
        model["Q", "arima"]  <- tmp$Q
        if(estimate){
          -logLik(model)
        } else model
      } else {
        if(estimate){
          1e100
        } else model
      }
    }
    fit_arima <- optim(par = c(rep(0, 4), log(1)), fn = likfn, method = "BFGS",
                       model = model_arima)
    print(fit_arima$par)
    model_arima <- likfn(fit_arima$par, model_arima, FALSE)
    out   <- KFS(model_arima)
    Out   <- ts(data.frame(u, obs = y,
                           smooth = signal(out)$signal,
                           prd    = signal(out,filtered = TRUE)$signal
    ), start = 1969, end = 2007, frequency = 1)
    ts.plot(window(Out, start = 1969, end = 2007),
            col = c("yellow",gray(0.5), "black", "blue"),
            ylim = c(-20,10),
            ylab = "u(t), feet per second", lty = c(1,3,1,1), lwd = c(6,2,2,2)
    )
    title(main = expression("Arima, "*dot(x)[t]))
    legend("bottomright",
           legend = c("Target", "Observed", "One-step ahead predictions", "Smoothed estimates" ),
           lty = c(1,3,1,1),
           lwd = c(6,2,2,2),
           col = c("gold", gray(0.5), "blue","black"),
           bty = "n")
    browser()

  }
}

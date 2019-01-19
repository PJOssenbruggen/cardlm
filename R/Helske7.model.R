#' \code{Helske7.model} is a wrapper function ring road experiment
#'
#' @param usd standard deviation of speed in mph, a number.
# #' @param zsd standard deviation of measuring location, a number.
#' @usage uses \code{Helske4.model} is an extension of \code{Helske} models 1, 2 and 3
#' @param veh, a number
# #' @examples
# #' Helske7.model(usd)
#' @export
Helske7.model <- function(usd,zsd) {
  set.seed(123)
  start  <- 0
  end    <- 40
  tseq   <- seq(0,40,0.125)
  Q0     <- 0
  r      <- 100                                   # road radius (feet)
  u      <- 78                                    # uniform speed fps
  z      <- u*tseq                                # distance travelled (feet)
  u      <- rep(u,length(tseq))
  u.e    <- rnorm(length(tseq), 0, usd*5280/3600) # Orbital speed noise
  v.e    <- u + u.e                               # v.e = orbital speed with noise
  z.e    <- rep(0, length(tseq))                  # distance travelled, orbit circumference location from starting point (x = r, y  =  0)
  z.e[1] <- 0
  for(i in 2:length(tseq)) z.e[i] <- z.e[i-1] + v.e[i]*0.125
  df     <- data.frame(t   = tseq, u, z, v.e, z.e)

  par(mfrow = c(2,2), pty = "s")
  plot(df$t, df$u, typ = "l", lwd = 2, col = "gold")
  lines(df$t, df$v.e, lty = 3, lwd = 1)
  title(main = "Speed")

  plot(df$t, df$z, typ = "l", lwd = 2, col = "gold", xlim = c(20,21), ylim = c(1540,1660))
  lines(df$t, df$z.e, lty = 3, lwd = 2)
  title(main = "Location")
  legend("topleft",
         legend = c("Gold Standard", "Observed"),
         lty = c(1,3),
         lwd = c(2,2),
         col = c("gold", gray(0.5)),
         bty = "n")
  browser()
  ########################################################################################
  alpha  <- sqrt((u^2-usd^2)/u^2)
  alpha  <- 1
  data   <- ts(df[,-1], start, end, frequency = 8)
  model  <- SSModel(data[,c(3,4)] ~  -1 +
                      SSMcustom(Z  = matrix(c(alpha,0.125,0,1),2,2),
                                T  = matrix(c(1,0.125,0,1),2,2),
                                R  = matrix(c(1,0.125/2),2,1),
                                Q  = matrix(NA),
                                P1inf = diag(2)
                      ),
                    H = matrix(c(NA,NA),2,2),
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

  update_model <- function(pars, model) {
    model["H"] <- pars[1:4]
    model["Q"] <- pars[5]
    model
  }
  check_model  <- function(model) (  model$H[1,1,1] > 0 &
                                     model$H[2,2,1] > 0 &
                                     model$Q[1,1,1] > 0
                                    )
  fit         <- fitSSM(model,
                        inits    = rep(0.1,17),
                        checkfn  = check_model,
                        updatefn = update_model,
                        method   = "BFGS")
  out         <- KFS(fit$model)
  print(out$P[,,end])

  Out <- ts(data.frame(df[,-1],
                       smooth.v = out$muhat[,1],
                       smooth.z = out$muhat[,2],
                       pred.v   = out$att[,1],
                       pred.z   = out$att[,2]
                       ),
            start, end, frequency = 8)

  ts.plot(window(Out, start = c(20,1), end = c(21,0))[,c(1,3,5,7)],
          col = c("gold", gray(0.5), "black", "blue"),
          ylab = "u(t), feet per second", lty = c(1,3,1,1), lwd = c(6,2,2,2)
  )
  title(main = "Speed Predictions")

  ts.plot(window(Out, start = c(20,1), end = c(21,0))[,c(2,4,6,8)],
          col = c("gold", gray(0.5), "black", "blue"),
          ylab = "x(t), feet", lty = c(1,3,1,1), lwd = c(6,2,2,2)
  )
  title(main = "Location Predictions")
  legend("topleft",
         legend = c("One-step ahead predictions", "Smoothed estimates" ),
         lty = c(1,1),
         lwd = c(2,2),
         col = c("blue","black"),
         bty = "n")
  browser()
  return(df)
}
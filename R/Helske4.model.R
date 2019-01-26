#' \code{Helske4.model} is a wrapper function straight road experiment
#'
#' @param usd standard deviation of speed in mph, a number.
#' @param zsd standard deviation of measuring location, a number.
#' @usage uses \code{Helske4.model} is an extension of \code{Helske} models 1, 2 and 3
#' @param veh, a number
# #' @examples
# #' Helske4.model(usd, zsd)
#' @export
Helske4.model <- function(usd,zsd) {
  set.seed(123)
  start  <- 0
  end    <- 40
  tseq   <- seq(0,40,0.125)
  Q0     <- 0
  u.e    <- rnorm(length(tseq), 0, usd*5280/3600) # Orbital speed noise
  u      <- 78            # uniform speed fps
  v.o    <- rep(u,length(tseq)) + u.e             # v.o = u + u.e, orbital speed with noise
  z.e    <- rnorm(length(tseq), 0, zsd)           # Deviation from the centerline (noise)
  r.o    <- rep(u,length(tseq)) + u.e             # v.o = u + u.e, orbital speed with noise
  r      <- 100                                   # road radius (feet)
  d      <- u*tseq                                # distance travelled, orbit circumference location from starting point (x = r, y  =  0)
  d.o    <- rep(0, length(tseq))
  for(i in 2:length(tseq)) d.o[i] <- d.o[i-1] + v.o[i]*0.125
  z.o    <- rep(r,length(tseq)) + z.e
  for(i in 2:length(tseq)) z.o[i] <- z.o[i-1] + z.e[i]
  df     <- data.frame(t = tseq, x = r*cos(d/r), y = r*sin(d/r),
                       u.x = rep(u,length(tseq))*cos(d/r), u.y = rep(u,length(tseq))*sin(d/r),
                       v.x = v.o*cos(d/r), v.y = v.o*sin(d/r),
                       z.x = z.o*cos(d.o/r), z.y = z.o*sin(d.o/r), d.o)
  df     <- cbind(df, theta.rad = asin(df$y/r))
##########################################################################################
  par(mfrow = c(2,2), pty = "s")
  plot(df$x,df$y, typ = "l", xlim = c(-120,120),
       ylim = c(-120,120), axes = FALSE, xlab = "", ylab = "", col = gray(0.8), lwd = 20)
  lines(df$x,df$y, col = "yellow")
  lines(df$z.x, df$z.y, cex = 0.25, col = "black", lty = 3)
  title(main = "Ring Road")
  arrows(0,0,100,0, angle = 25, length = 0.1)
  text(50,0, labels = "r = 50", pos = 3)

  data <- ts(df[,-1], start, end, frequency = 8)
  ts.plot(window(data, start = c(0,1), end = c(15,0))[,c(5,6)],
          col = c("black", "blue"),
          ylab = "z(t), feet", lty = c(1,1), lwd = c(2,2)
  )
  title(main = expression("("*dot(x)[t]*","*dot(y)[t]*")"))
  abline(h = 0, col = gray(0.5))

########################################################################################
  model  <- SSModel(data[,5] ~  -1 + SSMcycle(period = 2*pi, Q = NA),
                    H = matrix(NA,1,1),
                    distribution = "gaussian",
                    data = data,
                    tol  = .Machine$double.eps^0.5)
  update_model <- function(pars, model) {
    model["H"] <- pars[1]
    model["Q"] <- pars[2:3]
    model
  }
  check_model  <- function(model) (model["H"] > 0 &
                                     model$Q[1,1,1] > 0 &
                                     model$Q[2,2,1] > 0)
  fit         <- fitSSM(model,
                        inits    = rep(0.1,3),
                        checkfn  = check_model,
                        updatefn = update_model,
                        method   = "BFGS")
  out         <- KFS(fit$model)
  print(out$P[,,end])
  Out <- ts(data.frame(df[,-1], smooth = out$muhat, prd = out$att),start, end, frequency = 8)


  ts.plot(window(Out, start = c(0,1), end = c(20,0))[,c(3,5,11,12)],
          col = c("gold", gray(0.5), "black", "blue"),
          ylab = "u(t), feet per second", lty = c(1,3,1,1), lwd = c(6,2,2,2)
  )
  title(main = "SSMcycle Speed Predictions")



  p1 <- as.numeric(out$P[1,,][1,])
  p2 <- as.numeric(out$P[1,,][2,])
  plot(tseq, p1[-1], typ = "l", lwd = 2, ylim = c(0,max(c(p1,p2))))
  lines(tseq,p2[-1], lwd = 2)
  title(main = "Covariance")

########################################################################################
  if(FALSE) {
    model  <- SSModel(data[,c(3,4,5,6)] ~  -1 +
                        SSMcustom(Z  = matrix(c(1,0.125,0,0,0,1,0,0,0,0,1,0.125,0,0,0,1),4,4),
                                  T  = matrix(c(1,0.125,0,0,0,1,0,0,0,0,1,0.125,0,0,0,1),4,4),
                                  R  = matrix(c(1,0.125/2,1,0.125/2),4,1),
                                  Q  = matrix(NA),
                                  P1inf = diag(4)
                        ),
                      H = matrix(c(NA,NA,NA,NA),4,4),
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
    check_model  <- function(model) (model$H[1,1,1]  > 0 &
                                       model$H[2,2,1]  > 0 &
                                       model$H[3,3,1]  > 0 &
                                       model$H[4,4,1]  > 0 &
                                       model["Q"] > 0)
    update_model <- function(pars, model) {
      model["H"] <- pars[1:16]
      model["Q"] <- pars[17]
      model
    }
    fit         <- fitSSM(model,
                          updatefn = update_model,
                          checkfn  = check_model,
                          inits    = rep(0.1,17),
                          method   = "BFGS")
    out         <- KFS(fit$model)
    print(out$P[,,end])
    Out <- data.frame(df, smooth = out$muhat, prd = out$att)
    title(main = "Location/Speed Trend Model")
    lines(Out$smooth.z.x, Out$smooth.z.y, cex = 0.25, col = "black", lty = 1)
    lines(Out$prd.custom3, Out$prd.custom4, cex = 0.25, col = "blue", lty = 1)
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
    P  <- ts(out$P[1,,])
    ts.plot(window(P, start = c(0,1), end = c(15,0)), col = "black", ylab = "P", lwd = 2)
    title(main = "Covariance")

    ### Notes
    # 1. u = 50 = speed
    # 2. Q = 0 Here, the goal is to reach u = 50
    # 3. tracks well for all usd
    # 4. Covariance quickly reach steady state

    return(list(model, fit, out, df))
  }
}

#' \code{Helske6.model} is a wrapper function ring road experiment
#'
#' @param usd standard deviation of speed in mph, a number.
#' @param zsd standard deviation of measuring location, a number.
#' @usage uses \code{Helske4.model} is an extension of \code{Helske} models 1, 2 and 3
#' @param veh, a number
# #' @examples
# #' Helske6.model(usd,zsd)
#' @export
Helske6.model <- function(usd,zsd) {
  set.seed(123)
  start  <- 0
  end    <- 40
  tseq   <- seq(0,40,0.125)
  Q0     <- 0
  r      <- 100                                   # road radius (feet)
  d      <- 2*pi*r                                # road circumference
  u      <- 78                                    # uniform speed fps
  z      <- u*tseq                                # distance travelled (feet)
  u      <- rep(u,length(tseq))
  u.e    <- rnorm(length(tseq), 0, usd*5280/3600) # Orbital speed noise
  v.e    <- u + u.e                               # v.e = orbital speed with noise
  z.e    <- rep(0, length(tseq))                  # distance travelled, orbit circumference location from starting point (x = r, y  =  0)
  z.e[1] <- 0
  for(i in 2:length(tseq)) z.e[i] <- z.e[i-1] + v.e[i]*0.125
  theta  <- 2*z/r
  theta.e<- 2*z.e/r
  x      <- z*cos(theta)
  y      <- z*sin(theta)
  x.e    <- z*cos(theta.e)
  y.e    <- z*sin(theta.e)
  v.x    <- v.e*sin(theta.e)
  v.y    <- v.e*cos(theta.e)
  df     <- data.frame(t   = tseq, theta, theta.e, u, z, x, y, v.x, v.y, x.e, y.e,
                       degree   = 180/pi*theta, cosign = sign(cos(theta)), sinsign = sign(sin(theta)),
                       degree.e = 180/pi*theta.e, cosign.e = sign(cos(theta.e)), sinsign.e = sign(sin(theta.e)))

  df.check <- data.frame(u = df$u, u.check = sqrt(df$v.x^2 + df$v.y^2),
                         z = df$z, z.check = sqrt(df$x.e^2 + df$y.e^2))
  # u, z     = speed distance travelled along the ring road circumference with no noise
  # v.e, z.e = speed distance travelled along the ring road circumference with
#  print(head(df))
#  print(tail(df))

  par(mfrow = c(2,2), pty = "s")
  plot(r*cos(theta), r*sin(theta), typ = "l", xlim = c(-120,120),
       ylim = c(-120,120), axes = FALSE, xlab = "", ylab = "", col = gray(0.8), lwd = 20)
  lines(r*cos(theta), r*sin(theta), col = "white")
  for(i in 1:30) {
    points(r*cos(theta)[i], r*sin(theta)[i], cex = 0.75, col = "orange")
    points(r*cos(theta.e)[i], r*sin(theta.e)[i], cex = 0.5, col = "black", pch = 16)
    }
  title(main = "Ring Road")
  arrows(0,0,100,0, angle = 25, length = 0.1)
  text(50,0, labels = "r = 50", pos = 3)

  ##########################################################################################
  data <- ts(df[,-c(1:3)], start, end, frequency = 8)
  # print(df.check)
  if(FALSE) {
    ts.plot(window(data, start = c(0,1), end = c(40,0))[,c(5,6)],
            col = c("black", "blue"),
            ylab = "u(t), fps", lty = c(1,1), lwd = c(2,2)
    )
    abline(h = 0, col = gray(0.3))
    abline(v = 0, col = gray(0.3))
    title(main = expression("("*dot(x)[t]*","*dot(y)[t]*")"))
    ts.plot(window(data, start = c(0,1), end = c(40,0))[,c(7,8)],
            col = c("black", "blue"),
            ylab = "z(t), feet", lty = c(1,1), lwd = c(2,2)
    )
    abline(h = 0, col = gray(0.3))
    abline(v = 0, col = gray(0.3))
    title(main = expression("("*x[t]*","*y[t]*")"))

    ts.plot(window(sqrt(data[,7]^2+data[,8]^2), start = c(0,1), end = c(40,0)),
            col = c("black"),
            ylab = "z(t), feet", lty = c(1), lwd = c(2)
    )
    abline(h = 0, col = gray(0.3))
    abline(v = 0, col = gray(0.3))
    title(main = expression("("*x[t]*","*y[t]*") check"))
  }

  ########################################################################################
#  print(head(data))
#  print(tail(data))
  browser()
  alpha  <- sqrt((u^2-usd^2)/u^2)
  alpha  <- 1
  model  <- SSModel(data[,c(5,7,6,8)] ~  -1 +
                      SSMcustom(Z  = matrix(c(alpha,0.125,0,0,0,1,0,0,0,0,alpha,0.125,0,0,0,1),4,4),
                                T  = matrix(c(1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1),4,4),
                                R  = matrix(c(1,0.125/2,1,0.125/2),4,4),
                                Q  = matrix(NA,4,4),
                                P1inf = diag(4)
                      ),
                    H = matrix(NA,4,4),
                    distribution = "gaussian",
                    data = data,
                    tol  = .Machine$double.eps^0.5)
  update_model <- function(pars, model) {
    model["H"] <- pars[1:16]
    Q          <- diag(exp(pars[17:20]))
    Q[upper.tri(Q)] <- pars[21:26]
    model["Q"] <- crossprod(Q)
    model
  }
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
  check_model  <- function(model) (model$H[1,1,1] > 0 &
                                   model$H[2,2,1] > 0 &
                                   model$H[3,3,1] > 0 &
                                   model$H[4,4,1] > 0 &
                                   model$Q[1,1,1] > 0 &
                                   model$Q[2,2,1] > 0 &
                                   model$Q[3,3,1] > 0 &
                                   model$Q[4,4,1] > 0
                                   )
  fit         <- fitSSM(model,
                        inits    = rep(0.1,32),
                        checkfn  = check_model,
                        updatefn = update_model,
                        method   = "BFGS")
  out         <- KFS(fit$model, transform = "augment")
  print(out$P[,,end])
  print(out$Pinf)
browser()
  Out <- ts(data.frame(gold.v   = u,
                       gold.z   = z,
                       obs.v    = sqrt(data[,5]^2 + data[,6]^2),
                       obs.z    = sqrt(data[,7]^2 + data[,8]^2),
                       smooth.v = sqrt(out$muhat[,1]^2 + out$muhat[,3]^2),
                       smooth.z = sqrt(out$muhat[,2]^2 + out$muhat[,4]^2),
                       prd.v    = sqrt((out$att[,1] + out$att[,5])^2 +
                                         (out$att[,3] + out$att[,7])^2),
                       prd.x    = sqrt((out$att[,2] + out$att[,6])^2 +
                                         (out$att[,4] + out$att[,8])^2)
                       ),
                       start, end, frequency = 8)
  head(Out)
  tail(Out)
  browser()

  ts.plot(window(Out, start = c(0,1), end = c(40,1))[,c(1,3,5,7)],
          col = c("gold", gray(0.5), "black", "blue"),
          ylab = "u(t), feet per second", lty = c(1,3,1,1), lwd = c(6,2,2,2)
  )
  title(main = "SSMcycle Speed Model")

  ts.plot(window(Out, start = c(0,1), end = c(40,1))[,c(2,4,6,8)],
          col = c("gold", gray(0.5), "black", "blue"),
          ylab = "x(t), feet", lty = c(1,3,1,1), lwd = c(6,2,2,2)
  )
  title(main = "SSMcycle Location Model")

    p1 <- as.numeric(out$P[1,,][1,])
    p2 <- as.numeric(out$P[1,,][2,])
    plot(tseq, p1[-1], typ = "l", lwd = 2, ylim = c(0,max(c(p1,p2))), xlab = "Time", ylab = "P")
    lines(tseq,p2[-1], lwd = 2)
    title(main = "Covariance")
  browser()

  return(data)
  ########################################################################################
}

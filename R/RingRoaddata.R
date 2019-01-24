#' \code{RingRoaddata} is a wrapper function used for testing \code{Helske} models.
#'
#' @param usd standard deviation of speed in mph, a number.
#' @param zsd standard deviation of measuring location, a number.
# #' @examples
# #' RingRoaddata(usd)
#' @export
RingRoaddata <- function(usd,zsd) {
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
  df     <- data.frame(t   = tseq, theta, theta.e, u, z, v.e, z.e, x, y, v.x, v.y, x.e, y.e,
                       degree   = 180/pi*theta, cosign = sign(cos(theta)), sinsign = sign(sin(theta)),
                       degree.e = 180/pi*theta.e, cosign.e = sign(cos(theta.e)), sinsign.e = sign(sin(theta.e)))

  df.check <- data.frame(u = df$u, u.check = sqrt(df$v.x^2 + df$v.y^2),
                         z = df$z, z.check = sqrt(df$x.e^2 + df$y.e^2))

  df = df[,c(4:13,1:3,14)]

  layout(mat = matrix(c(1,1,2,3),2,2), widths = c(3,3), height = c(3,3))
  par(mar = c(1,3,1,3), pty = "s")
  plot(r*cos(theta), r*sin(theta), typ = "l", xlim = c(-120,120),
       ylim = c(-120,120), axes = FALSE, xlab = "", ylab = "", col = gray(0.8), lwd = 20)
  lines(r*cos(theta), r*sin(theta), col = "white")
  for(i in 1:30) {
    points(r*cos(theta)[i], r*sin(theta)[i], cex = 0.75, col = "gold", pch = 16)
    points(r*cos(theta.e)[i], r*sin(theta.e)[i], cex = 0.5, col = "black", pch = 16)
  }
  title(main = "Ring Road")
  arrows(0,0,100,0, angle = 25, length = 0.1)
  text(50,0, labels = "r = 50", pos = 3)
  legend(x = 0, y = -35, legend = c("target","observed"),
         col = c("gold","black"),
         pch = c(16,16), bty = "n")

  plot(tseq, u, col = "gold", ylab = "u(t), feet per second", lty = 1, lwd = 6,  xlab = "t, seconds")
  lines(tseq,v.e, col = "black", lty = 3, lwd = 2)
  title(main = expression(dot(x)[t]))
  legend("bottom", legend = c("target","observed"),
         col = c("gold","black"),
         lwd = c(6,2), lty = c(1,3), bty = "n")

  plot(tseq, z, col = "gold", ylab = "x(t), feet", lty = 1, lwd = 6,  xlab = "t, seconds")
  lines(tseq,z.e, col = "black", lty = 3, lwd = 2)
  title(main = expression(x[t]))

  return(df)
}

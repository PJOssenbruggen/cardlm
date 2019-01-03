#' \code{merge.analysis} is a wrapper function for estimating a successful, safe merge.
#'
#' @usage uses \code{merge.analysis}, \code{xab} and \code{uab} from the \code{cartools} Package.
# #' @examples
# #'
#' @export
merge.analysis <- function() {
  plot(tdf1df2[,1], tdf1df2[,3], typ = "l", ylim = c(-1200,2400), ylab = expression("x, feet"))
  xveh  <- seq(3,30,3)
  uveh  <- seq(2,29,3)
  tveh  <- tdf1df2[,1]
  h     <- x.xe <- u.x0  <- u.xe <- t.x0 <- t.xe <- rep(NA,10)
  u     <- tdf1df2[tdf1df2[,xveh[1]] <= 0, uveh[1]]
  t     <- tveh[tdf1df2[,xveh[1]] <= 0]
  u.x0[1] <- u[length(t)]
  t.x0[1] <- t[length(t)]
  x     <- tdf1df2[tdf1df2[,xveh[1]] <= -500, xveh[1]]
  u     <- tdf1df2[tdf1df2[,xveh[1]] <= -500, uveh[1]]
  t     <- tveh[tdf1df2[,xveh[1]] <= -500]
  x.xe[1] <- x[length(t)]
  u.xe[1] <- u[length(t)]
  t.xe[1] <- t[length(t)]
  for(i in 2:10) {
    lines(tdf1df2[,1], tdf1df2[,xveh[i]])
    u     <- tdf1df2[tdf1df2[,xveh[i]] <= 0, uveh[i]]
    t     <- tveh[tdf1df2[,xveh[i]] <= 0]
    u.x0[i] <- u[length(t)]
    t.x0[i] <- t[length(t)]
    x     <- tdf1df2[tdf1df2[,xveh[i]] <= -500, xveh[i]]
    u     <- tdf1df2[tdf1df2[,xveh[i]] <= -500, uveh[i]]
    t     <- tveh[tdf1df2[,xveh[i]] <= -500]
    x.xe[i] <- x[length(t)]
    u.xe[i] <- u[length(t)]
    t.xe[i] <- t[length(t)]
  }
  for(i in 1:10) h[i] <- hsafe(u.x0[i], 14)
  abline(h = c(0,-500))
  df <- data.frame(t.xe, u.xe, t.x0, u.x0, h, x.xe)
  for(i in 1:10) points(t.x0[i], -h[i])
  dt.0 <- t.x0 - t.xe
  df <- cbind(df, dt.0)
  for(i in 2:10) {
    lines(c(t.xe[i],t.x0[i]), c(x.xe[i], x.xe[i] + u.xe[i] * df[i,7]), lwd = 3, col = "red")
  }
  browser()
}

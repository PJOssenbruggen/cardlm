#' \code{mergedlm} is a wrapper function for plotting a ideal \code{zipper} and \code{side-by-side} bottleneck merges.
#'
#' @usage uses \code{xabparam}, \code{xab} and \code{uab} from the \code{cartools} Package.
# #' @examples
# #' mergedlm()
#' @export
mergedlm <- function() {
  tend  <- 12
  x1    <- 0
  u1    <- 53.1*5820/3600
  h     <- hsafe(u1,14)
  t1    <- 200/u1
  t2    <- 700/u1
  tseq  <- seq(0,tend,0.125)
  x1seq <- function(t) -700 + u1 * t
  x1seq <- x1seq(tseq)
  ab    <- xabparam(t1,t2,u1,u1,-500,-h)
  u2seq <- uab(u1, ab[1], ab[2], t1, t2)
  x2seq <- uab(u1, ab[1], ab[2], t1, t2)
  par(mfrow = c(1,2), pty = "s")
  plot(tseq,x1seq, xlab = "t, seconds", ylab = "x, feet", typ = "l",
       xlim = c(0,tend+0.5), ylim = c(-800,400), lwd = 2, col = "blue")
  abline(h = c(0,-500), col = gray(0.8))
  abline(v = 0, col = gray(0.8))
  abline(v = t2, col = gray(0.5), lty = 2)
  abline(v = t1, col = gray(0.5), lty = 2)
  points(t1,-500, pch = 16, cex = 1)
  points(t2,-h, pch = 16, cex = 1)
  veh2df <- data.frame(t = tseq, u = rep(0, length(tseq)), x = rep(0, length(tseq)))
  u2eff  <- (-h + 500)/(t2-t1)
  print(data.frame(u1,u2eff))
  tabseq <- tseq[tseq > t1 & tseq < t2]
  for(i in 1:length(tabseq)) veh2df[tabseq[i] == veh2df[,1],2] <- uab(u1,ab[1],ab[2], tabseq[i], t1)
  for(i in 1:length(tabseq)) veh2df[tabseq[i] == veh2df[,1],3] <- xab(-500,u1,ab[1],ab[2], tabseq[i], t1)
  x2 <- function(x0,u1,t, t1) x0 + u1* (t - t1)
  tabseq <- tseq[tseq <= t1]
  for(i in 1:length(tabseq)) veh2df[tabseq[i] == veh2df[,1],2] <- u1
  for(i in 1:length(tabseq)) veh2df[tabseq[i] == veh2df[,1],3] <- x2(-700,u1,tabseq[i],0)
  tabseq <- tseq[tseq >= t2]
  for(i in 1:length(tabseq)) veh2df[tabseq[i] == veh2df[,1],2] <- u1
  for(i in 1:length(tabseq)) veh2df[tabseq[i] == veh2df[,1],3] <- x2(-h,u1,tabseq[i],t2)
  lines(veh2df[,1],veh2df[,3], col = "red", lwd = 2)
  title(main = "Side-by-side Merge", sub = "A safe, 'ideal' merge.")
  axis(side = 3, at = t2, labels = expression(t[2]), line = -1, tick = FALSE)
  axis(side = 3, at = t1, labels = expression(t[1]), line = -1, tick = FALSE)
  axis(side = 4, at = 0, labels = expression(x[0]),  tick = TRUE)
  axis(side = 4, at = -500, labels = expression(x[e]), tick = TRUE)
  text(tend,x1seq[length(x1seq)], labels = 1, pos = 4, cex = 0.75)
  text(tend,veh2df[dim(veh2df)[1],3], labels = 2, pos = 4, cex = 0.75)
  plot(tseq,x1seq, xlab = "t, seconds", ylab = "x, feet", typ = "l",
       xlim = c(0,tend+0.5), ylim = c(-800,400), lwd = 2, col = "blue")
  abline(h = c(0,-500), col = gray(0.8))
  abline(v = 0, col = gray(0.8))
  abline(v = t2, col = gray(0.5), lty = 2)
  t1    <- (200+h/2)/u1
  points(t1,-500, pch = 16, cex = 1)
  points(t2,-h, pch = 16, cex = 1)
  abline(v = t1, col = gray(0.5), lty = 2)
  ab    <- xabparam(t1,t2,u1,u1,-500,-h)
  u2seq <- uab(u1, ab[1], ab[2], t1, t2)
  x2seq <- uab(u1, ab[1], ab[2], t1, t2)
  veh2df <- data.frame(t = tseq, u = rep(0, length(tseq)), x = rep(0, length(tseq)))
  u2eff  <- (-h + 500)/(t2-t1)
  print(data.frame(u1,u2eff))
  tabseq <- tseq[tseq > t1 & tseq < t2]
  for(i in 1:length(tabseq)) veh2df[tabseq[i] == veh2df[,1],2] <- uab(u1,ab[1],ab[2], tabseq[i], t1)
  for(i in 1:length(tabseq)) veh2df[tabseq[i] == veh2df[,1],3] <- xab(-500,u1,ab[1],ab[2], tabseq[i], t1)
  x2 <- function(x0,u1,t, t1) x0 + u1* (t - t1)
  tabseq <- tseq[tseq <= t1]
  for(i in 1:length(tabseq)) veh2df[tabseq[i] == veh2df[,1],2] <- u1
  for(i in 1:length(tabseq)) veh2df[tabseq[i] == veh2df[,1],3] <- x2(-700,u1,tabseq[i],0)
  tabseq <- tseq[tseq >= t2]
  x2 <- function(x0,u1,t, t1) x0 + u1* (t - t1)
  tabseq <- tseq[tseq <= t1]
  for(i in 1:length(tabseq)) veh2df[tabseq[i] == veh2df[,1],2] <- u1
  for(i in 1:length(tabseq)) veh2df[tabseq[i] == veh2df[,1],3] <- x2(-700-h/2,u1,tabseq[i],0)
  tabseq <- tseq[tseq >= t2]
  for(i in 1:length(tabseq)) veh2df[tabseq[i] == veh2df[,1],2] <- u1
  for(i in 1:length(tabseq)) veh2df[tabseq[i] == veh2df[,1],3] <- x2(-h,u1,tabseq[i],t2)
  lines(veh2df[,1],veh2df[,3], col = "red", lwd = 2)

  title(main = "Zipper Merge", sub = "A safe, 'ideal' merge.")
  axis(side = 3, at = t2, labels = expression(t[2]), line = -1, tick = FALSE)
  axis(side = 3, at = t1, labels = expression(t[1]), line = -1, tick = FALSE)
  axis(side = 4, at = 0, labels = expression(x[0]),  tick = TRUE)
  axis(side = 4, at = -500, labels = expression(x[e]), tick = TRUE)
  text(tend,x1seq[length(x1seq)], labels = 1, pos = 4, cex = 0.75)
  text(tend,veh2df[dim(veh2df)[1],3], labels = 2, pos = 4, cex = 0.75)
}

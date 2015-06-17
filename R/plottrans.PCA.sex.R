##' PCA plot for a transband, color by sex.
##'
##' @param sex Female =1, Male=2
##' @inheritParams plottrans.PCA
##' @export
plottrans.PCA.sex <- function(Y, sex, max.p=100, ...){

  ## use the first 100.
  Y <- Y[, 1:min(max.p, ncol(Y))]

  PC <- predict(princomp(Y))[, 1:2]

  Class <- sex
  xlim <- range(PC[, 1])
  ylim <- range(PC[, 2])
  px <- pretty(xlim)
  py <- pretty(ylim)

  blue <- rgb(123, 104, 238, maxColorValue = 256)
  orange <- rgb(230, 159, 0, maxColorValue = 256)
  green <- rgb(27, 159, 120, maxColorValue = 256)
  yellow <- rgb(255, 255, 0, maxColorValue = 256)
  color <- c(blue, orange, green, yellow)

  grayplot(x=PC[,1],y=PC[,2],
           pch=21,xat=px,yat=py,col="black",bg=color[Class],
           hlines=py, vlines=px,
           xaxt="n", yaxt="n",
           xaxs="r", yaxs="r",
           xlim=xlim, ylim=ylim,
           xlab="Principal Component 1", ylab="Principal Component 2",
           mgp=c(1.6,0.2,0), cex=0.8, ...)

  u <- par("usr")
  x <- u[1] + diff(u[1:2])*c(0.4, 0.6)
  points(x, rep(u[4]+diff(u[3:4])*0.035, 2), pch=21, col="black",
         bg=color, xpd=TRUE, cex=0.8)
  text(x, rep(u[4]+diff(u[3:4])*0.035, 2), c("Female","Male"),
       col=color, xpd=TRUE, cex=0.8, pos=4)
}

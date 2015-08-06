##' PCA plot for a transband.
##'
##' @param Y matrix of traits in the transband
##' @param geno genotype at the peak marker of scanone.mvn.
##' @param nonrecomb index for non-recombnent(in a 10cM region) individuals.
##' @param max.p max number of traits used for PCA analysis. This is
##' @param ... Optional graphics parameters passed to \code{\link[broman]{grayplot}}
##' used to avoid rank deficiency.
##' @export
plottrans.PCA <- function(Y, geno, nonrecomb, max.p=100, ...){

  ## use the first 100.
  Y <- Y[, 1:min(max.p, ncol(Y))]

  PC <- predict(princomp(Y))[, 1:2]

  blue <- rgb(123, 104, 238, maxColorValue = 256)
  orange <- rgb(230, 159, 0, maxColorValue = 256)
  green <- rgb(27, 159, 120, maxColorValue = 256)
  yellow <- rgb(255, 255, 0, maxColorValue = 256)
  genecolor <- c(blue, orange, green, yellow)

  Class <- geno
  Class[-nonrecomb] <- 4
  xlim <- range(PC[, 1])
  ylim <- range(PC[, 2])
  px <- pretty(xlim)
  py <- pretty(ylim)
  broman::grayplot(x=PC[,1],y=PC[,2],
                   pch=21,xat=px,yat=py,col="black",bg=genecolor[Class],
                   hlines=py, vlines=px,
                   xaxt="n", yaxt="n",
                   xaxs="r", yaxs="r",
                   xlim=xlim, ylim=ylim,
                   xlab="Principal Component 1", ylab="Principal Component 2",
                   mgp=c(1.6,0.2,0), cex=0.8, las=1, ...)

  u <- par("usr")
  x <- u[1] + diff(u[1:2])*((2:5)*0.1+0.05)
  points(x, rep(u[4]+diff(u[3:4])*0.035, 4), pch=21, col="black",
         bg=genecolor, xpd=TRUE, cex=0.8)
  x <- u[1] + diff(u[1:2])*((2:5)*0.1+0.05)
  text(x, rep(u[4]+diff(u[3:4])*0.035, 3), c("BB","BR","RR","Recombinants"),
       col=c(genecolor[1:3],"black"), xpd=TRUE, cex=0.8, pos=4)
}

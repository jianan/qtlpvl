##' LDA plot for a transband.
##'
##' This version of plot color the recombinant mice by their 2 QTL genotype instead of all yellow.
##'
##' @inheritParams plottrans.PCA
##' @export
plottrans.LDA2 <- function(Y, geno, nonrecomb, max.p=100, geno2, addlegend=TRUE, ...){
  suppressMessages(require(MASS))|| stop("the required package 'MASS' is not installed. ")

  ## use the first 100.
  Y <- Y[, 1:min(max.p, ncol(Y))]

  fit <- MASS::lda(Y[nonrecomb, ], grouping=geno[nonrecomb])
  pred <- predict(fit, Y[-nonrecomb, ])
  sca <- fit$scaling
  y.lda <- Y %*% sca

  orange <- rgb(230, 159, 0, maxColorValue = 256)
  green <- rgb(27, 159, 120, maxColorValue = 256)
  blue <- rgb(123, 104, 238, maxColorValue = 256)
  yellow <- rgb(255, 255, 0, maxColorValue = 256)
  genecolor <- c(blue,orange,green, yellow)

  Class <- geno
  Class[-nonrecomb] <- 4
  xlim <- range(y.lda[, 1])
  ylim <- range(y.lda[, 2])
  px <- pretty(xlim)
  py <- pretty(ylim)
  broman::grayplot(x=y.lda[,1],y=y.lda[,2],
                   pch=21,xat=px,yat=py,col="black",bg=genecolor[Class],
                   hlines=py, vlines=px,
                   xaxt="n", yaxt="n",
                   xaxs="r", yaxs="r",
                   xlim=xlim, ylim=ylim,
                   xlab="Linear Discriminant 1", ylab="Linear Discriminant 2",
                   mgp=c(1.6,0.2,0), cex=0.8, las=1, ...)
  if(addlegend){
    u <- par("usr")
    x <- u[1] + diff(u[1:2])*((2:5)*0.1+0.05)
    points(x, rep(u[4]+diff(u[3:4])*0.035, 4), pch=21, col="black",
           bg=genecolor, xpd=TRUE, cex=0.8)
    x <- u[1] + diff(u[1:2])*((2:5)*0.1+0.05)
    text(x, rep(u[4]+diff(u[3:4])*0.035, 3), c("BB","BR","RR","Recombinants"),
         col=c(genecolor[1:3],"black"), xpd=TRUE, cex=0.8, pos=4)
  }

  orange2 <- rgb(256, 159, 0, maxColorValue = 256)
  orange1 <- rgb(200, 159, 0, maxColorValue = 256)
  cols <- c(blue, "blue", "darkblue",
            orange1, orange, orange2,
            "green", "darkgreen", green)
  ## add cols for the recombinant by the two QTL model
  y.lda <- y.lda[-nonrecomb, ]
  geno2 <- geno2[-nonrecomb, ]
  col <- geno2[,1] + 3*(geno2[,2] - 1)
  if(addlegend){
    legends <- paste(rep(c("BB","BR","RR"),each=3),
                     rep(c("BB","BR","RR"),times=3),
                     sep=":")
    legend("bottomright", legend=legends, col=cols, pch=20, ncol=3)
  }
  points(y.lda, bg=cols[col], pch=21, cex=1, col="black")
  ## bring the informative points to the front
  ind <- which(!col %in% c(1,5,9))
  points(y.lda[ind,], bg=cols[col[ind]], pch=21, cex=1.3, col="black")
}

##' Plot LOD score versus QTL position for multiple traits.
##'
##' For each trait, plot its LOD score versus QTL position when the
##' LOD is bigger than a threshold.
##'
##' @inheritParams scanone.mvn
##' @param LOD.threshold threshold for QTL to be displayed.
##' @param ... Optional graphics arguments
##' @return A plot LOD score versus QTL position for multiple traits.
##'
##' @examples
##' set.seed(92950640)
##' data(listeria)
##' listeria <- calc.genoprob(listeria)
##' n <- nind(listeria)
##' chr <- "1"
##' geno <- pull.geno(listeria, chr=chr)
##' genotype1 <- geno[,7]
##' genotype2 <- geno[,10]
##' p <- 100
##' p1 <- floor(p/2)
##' G1 <- matrix(genotype1, n, p1)
##' G2 <- -matrix(genotype2, n, p-p1)
##' G2[G2==3] <- 2
##' G <- cbind(G1, G2*(-2))
##' Y <- matrix(rnorm(n*p),n,p)
##' Y <- Y + G
##' plotLOD(Y, listeria, chr)
##'
##' @export
plotLOD <- function(Y, cross, chr, addcovar=NULL, intcovar=NULL, LOD.threshold=3,  ...){

  n <- nrow(Y)
  p <- ncol(Y)
  p1 <- ncol(cross$pheno)
  cross$pheno <- data.frame(cross$pheno, Y)
  if(!is.null(colnames(Y))) names(cross$pheno)[p1+(1:p)] <- colnames(Y)
  out <- scanone(cross, pheno.col=p1+(1:p), method="hk", chr=chr,
                 addcovar=addcovar, intcovar=intcovar)
  maxPOS <- out$pos[apply(out[, -(1:2)], 2, which.max)]
  maxLOD <- apply(out[, -(1:2)], 2, max)

  x <- maxPOS
  y <- maxLOD
  isQTL <- maxLOD > LOD.threshold
  x <- x[isQTL]
  y <- y[isQTL]

  xlim <- range(x) * c(0.9, 1.1)
  ylim <- c(0, max(maxLOD)) * 1.1

  mgp <- c(1.6, 0.2, 0)
  plot(0, 0, type="n", mgp=mgp,
       xlim=xlim, ylim=ylim, ylab="", xlab="QTL pos (cM)",
       yaxt="n", xaxt="n", xaxs="i", las=1, tck=0, ...)
  bgcolor <- "gray80"
  u <- par("usr")
  rect(u[1], u[3], u[2], u[4], border=NA, col=bgcolor)
  abline(h=pretty(ylim, 10), v=pretty(xlim),col="white")
  axis(2, at=pretty(ylim,10),mgp=mgp,tick=FALSE,las=1)
  axis(1, at=pretty(xlim),mgp=mgp,tick=FALSE)
  title(ylab="LOD score", mgp=c(2.1, 0, 0))
  points(x=x, y=y, col=c("red", "blue")[ifelse(y>0, 1, 2)], pch=20, cex=0.7)
  abline(h=0)
  rect(u[1], u[3], u[2], u[4], border=TRUE)

  map <- unlist(pull.map(cross, chr=chr))
  map <- map[map>xlim[1] & map<xlim[2]]
  rug(map, ticksize=0.01)
}

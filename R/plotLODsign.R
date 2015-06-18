##' Plot the signed LOD score versus QTL position.
##'
##' For each trait, the QTL effect direction is used as the sign of
##' the LOD score.
##'
##' @inheritParams scanone.mvn
##' @param Y matrix, columns are quantitative traits maps to the same position.
##' @param LOD.threshold threshold for QTL to be displayed.
##' @param ... Optional graphics arguments
##' @return a plot the signed LOD score versus QTL position for multiple traits.
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
##' plotLODsign(Y, listeria, chr)
##'
##' @export
plotLODsign <- function(Y, cross, chr, LODsign, maxPOS, map,
                        addcovar=NULL, intcovar=NULL, LOD.threshold=3,
                        xlab="QTL pos (cM)", ylab="signed LOD score",
                        mgp=c(1.6, 0.2, 0), bgcolor="gray80",
                        ...){

  if(missing(LODsign) | missing(maxPOS)){
    n <- nrow(Y)
    p <- ncol(Y)
    p1 <- ncol(cross$pheno)
    cross$pheno <- data.frame(cross$pheno, Y)
    if(!is.null(colnames(Y))) names(cross$pheno)[p1+(1:p)] <- colnames(Y)
    if (!("prob" %in% names(cross[[c("geno",1)]]))){
      warning("First running calc.genoprob.")
      cross <- calc.genoprob(cross)
    }
    out <- scanone(cross, pheno.col=p1+(1:p), method="hk", chr=chr,
                   addcovar=addcovar, intcovar=intcovar)
    maxPOSind <- apply(out[, -(1:2)], 2, which.max)
    maxPOS <- out$pos[apply(out[, -(1:2)], 2, which.max)]
    maxLOD <- apply(out[, -(1:2)], 2, max)

    step <- attr(cross[[c("geno",chr,"prob")]],"step")
    off.end <- attr(cross[[c("geno",chr,"prob")]],"off.end")
    error.prob <- attr(cross[[c("geno",chr,"prob")]],"error.prob")
    map.function <- attr(cross[[c("geno",chr,"prob")]],"map.function")
    stepwidth <- attr(cross[[c("geno",chr,"prob")]],"stepwidth")
    geno <- argmax.geno(cross, step, off.end, error.prob, map.function, stepwidth)
    geno <- pull.argmaxgeno(geno, chr=chr)

    LODsign <- numeric(p)
    for(i in 1:p){
      LODsign[i] <- sign(mean(Y[geno[, maxPOSind[i]]==3, i]) - mean(Y[geno[, maxPOSind[i]]==1, i]))
    }
    LODsign <- maxLOD * LODsign
  }

  x <- maxPOS
  y <- LODsign
  isQTL <- abs(LODsign) > LOD.threshold
  x <- x[isQTL]
  y <- y[isQTL]

  m <- max(abs(LODsign))
  xlim <- range(x) * c(0.9, 1.1)
  ylim <- c(-m, m) * 1.1
  xp <- pretty(xlim, 10)
  yp <- pretty(ylim, 10)

  plot(0, 0, type="n", mgp=mgp,
       xlim=xlim, ylim=ylim, ylab="", xlab=xlab,
       yaxt="n", xaxt="n", xaxs="i", las=1, ...)
  u <- par("usr")
  rect(u[1], u[3], u[2], u[4], border="black", col=bgcolor)
  abline(h=yp, v=xp, col="white")
  axis(2, at=yp, mgp=mgp, tick=FALSE, las=1)
  axis(1, at=xp, mgp=mgp, tick=FALSE)
  title(ylab=ylab, mgp=c(2.1, 0, 0))
  points(x=x, y=y, col=c("violetred", "slateblue")[ifelse(y>0, 1, 2)], pch=20, cex=0.7)
  abline(h=0)
  box()

  if(missing(map) & !missing(cross) & !missing(chr)){
    map <- pull.map(cross, chr=chr)[[1]]
  }
  map <- map[map > xlim[1] & map < xlim[2]]
  rug(map, ticksize=0.01)

}

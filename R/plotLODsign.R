##' Plot the signed LOD score versus QTL position.
##'
##' For each trait, the QTL effect direction is used as the sign of
##' the LOD score. 
##' 
##' @inheritParams scanone.mvn
##' @param Y matrix, columns are quantitative traits maps to the same position. 
##' @param LOD.threshold threshold for QTL to be displayed.
##' @param ... Optional graphics arguments
##' @return a plot of QTL dominance effect versus QTL additive effect.
##' @export
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
##' plotLODsign(listeria, Y, chr)

plotLODsign <- function(cross, Y, chr, addcovar=NULL, intcovar=NULL, LOD.threshold=3,  ...){
  n <- nrow(Y)
  p <- ncol(Y)
  p1 <- ncol(cross$pheno)
  cross$pheno <- data.frame(cross$pheno, Y)
  if(!is.null(colnames(Y))) names(cross$pheno)[p1+(1:p)] <- colnames(Y)
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
    LODsign[i] <- sign(mean(Y[geno[, maxPOSind[i]]==1, i]) - mean(Y[geno[, maxPOSind[i]]==3, i]))
  }

  x <- maxPOS
  y <- maxLOD * LODsign
  isQTL <- maxLOD > LOD.threshold
  x <- x[isQTL]
  y <- y[isQTL]

  xlim <- range(x) * c(0.9, 1.1)
  ylim <- c(-max(maxLOD), max(maxLOD)) * 1.1
  
  bgcolor <- "gray80"
  grayplot(0, 0, type="n", mgp=c(1.6, 0.2, 0), 
           yat=pretty(ylim, 10), xat=pretty(xlim, 10), xlim=xlim, ylim=ylim, 
           hlines=pretty(ylim, 10), vlines=pretty(x), ylab="", xlab="QTL pos (cM)", 
           yaxt="n", xaxt="n", xaxs="i", bgcolor=bgcolor, ...)
  title(ylab="signed LOD score", mgp=c(2.1, 0, 0))
  points(x=x, y=y, col=c("red", "blue")[ifelse(y>0, 1, 2)], pch=20, cex=0.7)
  abline(h=0)
  
}



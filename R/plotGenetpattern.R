##' Plot the QTL dominance effect versus QTL additive effect.
##'
##' QTL additive effect is defined as a=(RR-BB)/2;
##' QTL dominant effect is defined as d=BR-(BB+RR)/2.
##'
##' @inheritParams scanone.mvn
##' @inheritParams plotLODsign
##' @param genotype QTL genotype for the common QTL.
##' @param a vector of additive effect
##' @param d vector of diminance effect
##' @param xlab x-axis label
##' @param ylab y-axis label
##' @param mgp Margin parameters
##' @return a plot of QTL dominance effect versus QTL additive effect.
##'
##' @examples
##' data(fake.phenos)
##' data(listeria)
##' listeria <- calc.genoprob(listeria, step=1)
##' plotGenetpattern(Y=fake.phenos, cross=listeria, chr=1)
##'
##' @export
plotGenetpattern <- function(Y, genotype=NULL, cross, chr, a, d,
                             addcovar=NULL, intcovar=NULL, LOD.threshold=3,
                             xlab="Additive effect",
                             ylab= "Dominance effect",
                             mgp=c(1.6, 0.2, 0),
                             ...){

  if(missing(a) | missing(d)){
    if(!is.matrix(Y)) stop("Y need to be a matrix of quantitative traits")
    if(!missing(genotype) & length(genotype) != nrow(Y))
        stop("number of objects in 'Y' and 'genotype' not same. ")
    if(!missing(genotype) & length(unique(genotype)) != 3)
        stop("genotype should have 3 levels. ")

    if(missing(genotype)){  ## if genotype is not given, use QTL for each trait.

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
      maxLOD <- apply(out[, -(1:2)], 2, max)

      step <- attr(cross[[c("geno",chr,"prob")]],"step")
      off.end <- attr(cross[[c("geno",chr,"prob")]],"off.end")
      error.prob <- attr(cross[[c("geno",chr,"prob")]],"error.prob")
      map.function <- attr(cross[[c("geno",chr,"prob")]],"map.function")
      stepwidth <- attr(cross[[c("geno",chr,"prob")]],"stepwidth")
      geno <- argmax.geno(cross, step, off.end, error.prob, map.function, stepwidth)
      geno <- pull.argmaxgeno(geno, chr=chr)

      eff1 <- eff2 <- eff3 <- numeric(p)
      for(i in 1:p){
        eff1[i] <- mean(Y[geno[, maxPOSind[i]]==1, i], na.rm=TRUE)
        eff2[i] <- mean(Y[geno[, maxPOSind[i]]==2, i], na.rm=TRUE)
        eff3[i] <- mean(Y[geno[, maxPOSind[i]]==3, i], na.rm=TRUE)
      }
      eff1[maxLOD <= LOD.threshold] <- NA
      eff2[maxLOD <= LOD.threshold] <- NA
      eff3[maxLOD <= LOD.threshold] <- NA
    }else{
      gn <- as.numeric(genotype)
      eff1 <- apply(Y[which(gn==1), ], 2, mean, na.rm=TRUE)
      eff2 <- apply(Y[which(gn==2), ], 2, mean, na.rm=TRUE)
      eff3 <- apply(Y[which(gn==3), ], 2, mean, na.rm=TRUE)
    }
    Eff <- cbind(eff1, eff2, eff3)
    a <- (Eff[, 3] - Eff[, 1])/2
    d <- Eff[, 2] - (Eff[, 3] + Eff[, 1])/2
  }

  lim <- max(abs(c(a,d)), na.rm=TRUE)
  if(lim == -Inf) lim <- 1
  xlim <- ylim <- c(-lim, lim)*1.1
  plot(x=a, y=d, pch=20, xlim=xlim, ylim=ylim, xlab="", ylab=ylab, type="n",
       mgp=mgp, las=1, tck=0, xaxt="n", yaxt="n", ...)
  axis(side=1, mgp=mgp+c(0,0.2,0), tick=TRUE, tck=-0.02)
  axis(side=2, mgp=mgp+c(0,0.2,0), tick=TRUE, tck=-0.02)
  title(xlab=xlab, mgp=mgp-c(0.5,0,0), ...)
  u <- par("usr")
  rect(u[1], u[3], u[2], u[4], col="gray80")
  abline(v = 0, h = 0, col="white")
  abline(a = 0, b = -1, col="white")
  abline(a = 0, b = +1, col="white")
  points(x=a, y=d, col=c("violetred", "slateblue")[ifelse(a>0, 1, 2)], pch=20, cex=0.7)
  rect(u[1], u[3], u[2], u[4])
}

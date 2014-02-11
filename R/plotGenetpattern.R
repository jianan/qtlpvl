##' Plot the QTL dominance effect versus QTL additive effect.
##'
##' QTL additive effect is defined as a=(RR-BB)/2;
##' QTL dominant effect is defined as d=BR-(BB+RR)/2.
##' 
##' @param Y matrix, columns are quantitative traits maps to the same position. 
##' @param genotype QTL genotype for the traits.
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
##' plotGenetpattern(Y, genotype1)

plotGenetpattern <- function(Y, genotype){
  if(!is.matrix(Y)) stop("Y need to be a matrix of quantitative traits")
  if(length(genotype) != nrow(Y)) stop("number of objects in 'Y' and 'genotype' not same. ")
  if(length(unique(genotype)) != 3) stop("genotype should have 3 levels. ")
  
  gn <- as.numeric(genotype)
  eff1 <- apply(Y[which(gn==1), ], 2, mean)
  eff2 <- apply(Y[which(gn==2), ], 2, mean)
  eff3 <- apply(Y[which(gn==3), ], 2, mean)
  Eff <- cbind(eff1, eff2, eff3)
  a <- (Eff[, 3] - Eff[, 1])/2
  d <- Eff[, 2] - (Eff[, 3] + Eff[, 1])/2 
  lim <- max(abs(c(a,d)))
  xlim <- ylim <- c(-lim, lim)*1.1
  plot(x=a, y=d, pch=20, xlim=xlim, ylim=ylim, 
       xlab="additive effect: a=(RR-BB)/2", ylab= "dominance effect: d=BR-(BB+RR)/2", 
       mgp=c(1.6, 0.2, 0))
  abline(v=0, h=0, a=0, b=1)
  abline(a=0, b=-1)

}


\name{plotGenetpattern}
\alias{plotGenetpattern}
\title{Plot the QTL dominance effect versus QTL additive effect.}
\usage{
plotGenetpattern(Y, genotype = NULL, cross, chr, addcovar = NULL,
  intcovar = NULL, LOD.threshold = 3, ...)
}
\arguments{
  \item{genotype}{QTL genotype for the common QTL.}

  \item{Y}{matrix of multiple traits with n rows.}

  \item{cross}{An object of class \code{cross}. See
  \code{read.cross} for details.}

  \item{chr}{Optional vector indicating the chromosomes for
  which LOD scores should be calculated.  This should be a
  vector of character strings referring to chromosomes by
  name; numeric values are converted to strings.}

  \item{addcovar}{Additive covariates.}

  \item{intcovar}{Interactive covariates.}

  \item{LOD.threshold}{threshold for QTL to be displayed.}

  \item{...}{Optional graphics arguments}
}
\value{
a plot of QTL dominance effect versus QTL additive effect.
}
\description{
QTL additive effect is defined as a=(RR-BB)/2; QTL dominant
effect is defined as d=BR-(BB+RR)/2.
}
\examples{
set.seed(92950640)
data(listeria)
listeria <- calc.genoprob(listeria)
n <- nind(listeria)
chr <- "1"
geno <- pull.geno(listeria, chr=chr)
genotype1 <- geno[,7]
genotype2 <- geno[,10]
p <- 100
p1 <- floor(p/2)
G1 <- matrix(genotype1, n, p1)
G2 <- -matrix(genotype2, n, p-p1)
G2[G2==3] <- 2
G <- cbind(G1, G2*(-2))
Y <- matrix(rnorm(n*p),n,p)
Y <- Y + G
plotGenetpattern(Y, genotype=genotype1)
plotGenetpattern(Y, cross=listeria, chr=chr)
}


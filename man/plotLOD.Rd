\name{plotLOD}
\alias{plotLOD}
\title{Plot LOD score versus QTL position for multiple traits.}
\usage{
plotLOD(cross, Y, chr, addcovar = NULL, intcovar = NULL,
  LOD.threshold = 3, ...)
}
\arguments{
  \item{LOD.threshold}{threshold for QTL to be displayed.}

  \item{...}{Optional graphics arguments}

  \item{cross}{An object of class \code{cross}. See
  \code{read.cross} for details.}

  \item{Y}{matrix of multiple traits with n rows.}

  \item{chr}{Optional vector indicating the chromosomes for
  which LOD scores should be calculated.  This should be a
  vector of character strings referring to chromosomes by
  name; numeric values are converted to strings.}

  \item{addcovar}{Additive covariates.}

  \item{intcovar}{Interactive covariates.}
}
\value{
A plot LOD score versus QTL position for multiple traits.
}
\description{
For each trait, plot its LOD score versus QTL position when
the LOD is bigger than a threshold.
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
plotLOD(listeria, Y, chr)
}


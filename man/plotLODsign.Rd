\name{plotLODsign}
\alias{plotLODsign}
\title{Plot the signed LOD score versus QTL position.}
\usage{
plotLODsign(cross, Y, chr, addcovar = NULL, intcovar = NULL,
  LOD.threshold = 3, ...)
}
\arguments{
  \item{Y}{matrix, columns are quantitative traits maps to
  the same position.}

  \item{LOD.threshold}{threshold for QTL to be displayed.}

  \item{...}{Optional graphics arguments}

  \item{cross}{An object of class \code{cross}. See
  \code{read.cross} for details.}

  \item{chr}{Optional vector indicating the chromosomes for
  which LOD scores should be calculated.  This should be a
  vector of character strings referring to chromosomes by
  name; numeric values are converted to strings.}

  \item{addcovar}{Additive covariates.}

  \item{intcovar}{Interactive covariates.}
}
\value{
a plot the signed LOD score versus QTL position for
multiple traits.
}
\description{
For each trait, the QTL effect direction is used as the
sign of the LOD score.
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
plotLODsign(listeria, Y, chr)
}


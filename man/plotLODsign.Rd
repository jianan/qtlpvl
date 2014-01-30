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

  \item{genotype}{QTL genotype for the traits.}
}
\value{
a plot of QTL dominance effect versus QTL additive effect.
}
\description{
For each trait, the QTL effect direction is used as the
sign of the LOD score.
}
\examples{
data(listeria)
listeria <- calc.genoprob(listeria)
geno <- pull.geno(listeria, chr="1")
genotype1 <- geno[,7]
genotype2 <- geno[,10]
n <- length(genotype1)
p <- 100
p1 <- floor(p/2)
G1 <- matrix(genotype1, n, p1)
G2 <- -matrix(genotype2, n, p-p1)
G <- cbind(G1, G2)
Y <- matrix(rnorm(n*p,sd=0.5),n,p)
Y <- Y + G
plotLODsign(listeria, Y, chr="1")
}


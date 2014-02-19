\name{scanone.mvn}
\alias{scanone.mvn}
\title{Multivariate normal genome scan with a single QTL model.}
\usage{
scanone.mvn(cross, Y, chr = NULL, addcovar = NULL, intcovar = NULL,
  tol = 1e-07)
}
\arguments{
  \item{cross}{An object of class \code{cross}. See
  \code{read.cross} for details.}

  \item{Y}{matrix of multiple traits with n rows.}

  \item{chr}{Optional vector indicating the chromosomes for
  which LOD scores should be calculated.  This should be a
  vector of character strings referring to chromosomes by
  name; numeric values are converted to strings.}

  \item{addcovar}{Additive covariates.}

  \item{intcovar}{Interactive covariates.}

  \item{tol}{Tolerance value for the \code{qr}
  decomposition in \code{lm} fitting.}
}
\value{
A data.frame whose first column contains the chromosome
IDs, second column contains cM positions, third column
contains the LOD scores.
}
\description{
Genome scan with a single QTL model, with allowance for
possible covariates, using multivariate normal model for
the multiple traits.
}
\examples{
set.seed(92950640)
data(listeria)
listeria <- calc.genoprob(listeria,step=1)
n <- nind(listeria)
chr <- "1"
geno <- pull.geno(listeria, chr=chr)
genotype1 <- geno[,7]
genotype2 <- geno[,10]
p <- 10
p1 <- floor(p/2)
G1 <- matrix(genotype1, n, p1)
G2 <- matrix(genotype2, n, p-p1)
G2[G2==3] <- 2
G <- cbind(G1, G2*(-2))
Y <- matrix(rnorm(n*p),n,p)
Y <- Y + G
out <- scanone.mvn(listeria, Y, chr)
summary(out)
plot(out)
}
\seealso{
qtl::scanone
}
\keyword{QTL}


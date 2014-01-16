\name{scanone.mvn}
\alias{scanone.mvn}
\title{Multivariate normal genome scan with a single QTL model}
\usage{
scanone.mvn(Y, cross, chr = NULL, addcov = NULL, intcov = NULL,
  tol = 1e-07)
}
\arguments{
  \item{Y}{matrix of multiple traits with n rows.}

  \item{cross}{An object of class \code{cross}. See
  \code{read.cross} for details.}

  \item{chr}{Optional vector indicating the chromosomes for
  which LOD scores should be calculated.  This should be a
  vector of character strings referring to chromosomes by
  name; numeric values are converted to strings.}

  \item{addcov}{Additive covariates.}

  \item{intcov}{Interactive covariates.}

  \item{tol}{Tolerance value for the \code{qr}
  decomposition in \code{lm} fitting.}
}
\value{
A data.frame whose first column contains the chromosome
IDs, second column contains cM positions, third column
contains the LOD scores.
}
\description{
Genome scan with a single QTL model, with possible
allowance for covariates, using multivariate normal model
for the multiple traits.
}
\examples{
library(qtl)
data(hyper)
n <- 250
p <- 5
Y <- matrix(rnorm(n*p),n,p)
scanone.mvn(Y, hyper)
}
\seealso{
qtl::scanone
}
\keyword{QTL}


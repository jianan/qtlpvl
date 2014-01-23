\name{scanone.mvn}
\alias{scanone.mvn}
\title{Multivariate normal genome scan with a single QTL model.}
\usage{
scanone.mvn(cross, Y, chr = NULL, addcovar = NULL, intcovar = NULL,
  tol = 0.0000001)
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
data(hyper)
n <- 250
p <- 5
Y <- matrix(rnorm(n*p),n,p)
hyper <- calc.genoprob(hyper)
scanone.mvn(hyper, Y)
}
\seealso{
qtl::scanone
}
\keyword{QTL}


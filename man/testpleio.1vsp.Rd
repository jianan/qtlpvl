\name{testpleio.1vsp}
\alias{testpleio.1vsp}
\title{Test of pleiotrophic QTL vs close linked QTLs.}
\usage{
testpleio.1vsp(cross, Y, chr = "6", addcov = NULL, intcov = NULL,
  n.simu = NA, tol = 1e-07)
}
\arguments{
  \item{Y}{matrix of multiple traits with rows equals to n,
  the number of individuals.}

  \item{cross}{An object of class \code{cross}. See
  \code{read.cross} for details.}

  \item{chr}{\code{chr, region.l, region.r} are used to
  specify the interval of interest. see also
  \code{int.method}.}

  \item{addcov}{Additive covariates.}

  \item{intcov}{Interactive covariates.}

  \item{n.simu}{number of simulations for p-value.}

  \item{tol}{Tolerance value for the \code{qr}
  decomposition in \code{lm} fitting.}
}
\value{
a list of LOD1, LODp, LODdiff, LOD1.pos and maxPOS ...
P-value
}
\description{
Test if multitriats are controlled by one pleiotrophic QTL
or multiple close linked QTLs, allowing each trait to have
its own QTL. p-value is calculated from paramatric
bootstrap, which takes the estimated paramaters under null
hypothesis and generate data from them to get the null
empirical distribution of test statistic.
}
\examples{
library(qtl)
data(hyper)
n <- 250
p <- 5
Y <- matrix(rnorm(n*p),n,p)
testpleio.1vsp(cross=hyper, Y=Y, chr="2")
}


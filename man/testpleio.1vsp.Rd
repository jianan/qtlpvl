\name{testpleio.1vsp}
\alias{testpleio.1vsp}
\title{Test of pleiotrophic QTL vs close linked QTLs.}
\usage{
testpleio.1vsp(cross, Y, chr = "6", addcovar = NULL, intcovar = NULL,
  n.simu = NA, tol = 1e-07)
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
data(hyper)
n <- 250
p <- 5
Y <- matrix(rnorm(n*p),n,p)
hyper <- calc.genoprob(hyper)
summary(testpleio.1vsp(cross=hyper, Y=Y, chr="2"))
}


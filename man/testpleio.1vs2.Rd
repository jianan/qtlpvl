\name{testpleio.1vs2}
\alias{testpleio.1vs2}
\title{Test of one pleiotrophic QTL vs two close linked QTLs.}
\usage{
testpleio.1vs2(cross, Y, chr = "6", addcov = NULL, intcov = NULL,
  region.l = NA, region.r = NA, int.method = "bayes", search = "fast",
  RandomStart = TRUE)
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

  \item{region.l}{left bound}

  \item{region.r}{right bound}

  \item{int.method}{"bayes" or "1.5lod" method to
  calculated the interval of interest if \code{region.l}
  and \code{region.r} is not specified.}

  \item{search}{searching method for two-QTL model.}

  \item{RandomStart}{use random starting point for the
  two-QTL model or not. default is \code{TRUE}.}
}
\value{
a list of sth
}
\description{
Test if multitriats are controlled by one pleiotrophic QTL
or two close linked QTLs. p-value is calculated from
paramatric bootstrap, which takes the estimated paramaters
under null hypothesis and generate data from them to get
the null empirical distribution of test statistic.
}
\examples{
library(qtl)
data(hyper)
n <- 250
p <- 5
Y <- matrix(rnorm(n*p),n,p)
testpleio.1vs2(cross=hyper, Y=Y, chr="2")
}


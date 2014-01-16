\name{test.pleio.vs.2qtl}
\alias{test.pleio.vs.2qtl}
\title{Test of one pleiotrophic QTL vs two close linked QTLs.}
\usage{
test.pleio.vs.2qtl(cross, Y, chr = "6", addcov = NULL, intcov = NULL,
  region.l = NA, region.r = NA, int.method = "bayes", search = "fast",
  RandomStart = TRUE)
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

  \item{region.l}{left bound}

  \item{region.r}{right bound}

  \item{int.method}{"bayes" or "1.5lod"}

  \item{search}{searching method}

  \item{RandomStart}{default is TRUE}
}
\value{
a list of sth
}
\description{
Test if multitriats are controlled by one pleiotrophic QTL
or two close linked QTLs. p-value is calculated from
paramatric bootstrap.
}
\examples{
library(qtl)
data(hyper)
n <- 250
p <- 5
Y <- matrix(rnorm(n*p),n,p)
scanone.mvn(Y, hyper)
}


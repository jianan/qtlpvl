\name{testpleio.1vs2}
\alias{testpleio.1vs2}
\title{Test of one pleiotrophic QTL vs two close linked QTLs.}
\usage{
testpleio.1vs2(cross, Y, chr = "6", addcovar = NULL, intcovar = NULL,
  region.l = NA, region.r = NA, int.method = "bayes", search = "fast",
  RandomStart = TRUE, n.simu = 1000, tol = 1e-07)
}
\arguments{
  \item{region.l}{left bound}

  \item{region.r}{right bound}

  \item{int.method}{"bayes" or "1.5lod" method to
  calculated the interval of interest if \code{region.l}
  and \code{region.r} is not specified.}

  \item{search}{searching method for two-QTL model.}

  \item{RandomStart}{use random starting point for the
  two-QTL model or not. default is \code{TRUE}.}

  \item{n.simu}{number of simulations for p-value.}

  \item{tol}{Tolerance value for the \code{qr}
  decomposition in \code{lm} fitting.}

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
a list of LODdiff, Group, P-value... \item{LOD1}{One
dimesional joint mapping.} \item{LOD2}{LOD score for
estimate of best two QTL model, each trait is influenced by
the left or the right QTL.} \item{LODdiff}{Difference of
best one QTL model and best two QTL model, LODdiff =
max(LOD2) - max(LOD1)} \item{pvalue}{P-value from
parametric bootstrap simulations.}
}
\description{
Test if multitriats are controlled by one pleiotrophic QTL
or two close linked QTLs. p-value is calculated from
paramatric bootstrap, which takes the estimated paramaters
under null hypothesis and generate data from them to get
the null empirical distribution of test statistic.
}
\examples{
data(hyper)
n <- 250
p <- 5
Y <- matrix(rnorm(n*p),n,p)
hyper <- calc.genoprob(hyper)
summary(testpleio.1vs2(cross=hyper, Y=Y, chr="2", n.simu=100))
}


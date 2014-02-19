\name{testpleio.1vs2}
\alias{testpleio.1vs2}
\title{Test of one pleiotrophic QTL vs two close linked QTLs.}
\usage{
testpleio.1vs2(cross, Y, chr = "6", addcovar = NULL, intcovar = NULL,
  region.l = NA, region.r = NA, int.method = "bayes", search = "fast",
  RandomStart = TRUE, RandomCut = FALSE, simu = "parametric",
  n.simu = 1000, tol = 1e-07)
}
\arguments{
  \item{region.l}{left bound}

  \item{region.r}{right bound}

  \item{int.method}{"bayes" or "1.5lod" method to
  calculated the interval of interest if \code{region.l}
  and \code{region.r} is not specified.}

  \item{search}{searching method for two-QTL model, "fast"
  or "complete".}

  \item{RandomStart}{use random starting point for the
  two-QTL model or not. default is \code{TRUE}.}

  \item{RandomCut}{use random cutting or not when there are
  traits mapped to the same location. default is "FALSE".}

  \item{simu}{"parametric" or "permutation" method for
  simulations. default is "parametric".}

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
a list. \item{LOD1}{LOD score of one dimensional joint
mapping.} \item{LOD2}{LOD score for estimate of best two
QTL model, each trait is influenced by the left or the
right QTL.} \item{LODdiff}{Difference of best one QTL model
and best two QTL model, LODdiff = max(LOD2) - max(LOD1)}
\item{pvalue}{P-value from parametric bootstrap
simulations.} \item{Group}{Indicating which traits are
influenced by which QTL under alternative.}
}
\description{
Test if multitriats are controlled by one pleiotrophic QTL
or two close linked QTLs. p-value is calculated from
parametric bootstrap, which takes the estimated parameters
under null hypothesis and generate data from them to get
the null empirical distribution of test statistic.
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
obj <- testpleio.1vs2(listeria, Y, chr, n.simu=100,
                      region.l=60, region.r=90)
summary(obj)
plot(obj)
}


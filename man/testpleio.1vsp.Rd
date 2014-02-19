\name{testpleio.1vsp}
\alias{testpleio.1vsp}
\title{Test of pleiotrophic QTL vs close linked QTLs.}
\usage{
testpleio.1vsp(cross, Y, chr = "6", addcovar = NULL, intcovar = NULL,
  n.simu = 1000, tol = 1e-07)
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
P-value \item{LOD1}{LOD score of one dimensional joint
mapping.} \item{LODp}{LOD score for estimate of p QTL
model, each trait is influenced by a separate QTL.}
\item{LODdiff}{Difference of best one QTL model and best p
QTL model, LODdiff = max(LODp) - max(LOD1)}
\item{pvalue}{P-value from parametric bootstrap
simulations.}
}
\description{
Test if multitriats are controlled by one pleiotrophic QTL
or multiple close linked QTLs, allowing each trait to have
its own QTL. p-value is calculated from parametric
bootstrap, which takes the estimated parameters under null
hypothesis and generate data from them to get the null
empirical distribution of test statistic.
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
obj <- testpleio.1vsp(listeria, Y, chr, n.simu=100)
summary(obj)
plot(obj)
}


##' Test of pleiotrophic QTL vs close linked QTLs.
##'
##' Test if multitriats are controlled by one pleiotrophic QTL or
##' multiple close linked QTLs, allowing each trait to have its own
##' QTL. p-value is calculated from paramatric bootstrap, which takes
##' the estimated paramaters under null hypothesis and generate data
##' from them to get the null empirical distribution of test
##' statistic.
##' 
##' @param Y matrix of multiple traits with rows equals to n, the
##' number of individuals.
##' @param cross An object of class \code{cross}. See
##' \code{read.cross} for details.
##' @param chr \code{chr, region.l, region.r} are used to specify the
##' interval of interest. see also \code{int.method}.
##' @param addcov Additive covariates.
##' @param intcov Interactive covariates.
##' @return a list of LOD1, LODp, LODdiff, LOD1.pos and maxPOS ...
##' 
##' @export
##' @examples
##' library(qtl)
##' data(hyper)
##' n <- 250
##' p <- 5
##' Y <- matrix(rnorm(n*p),n,p)
##' testpleio.1vsp(cross=hyper, Y=Y, chr="2")

testpleio.1vsp <- function(cross, Y, chr="6", addcov=NULL, intcov=NULL, tol=1e-7){

  ## - scanone for each trait to get E.matrix
  ## - 1 qtl model: same qtl for all trait
  ## - LODdiff = max(LODp) - max(LOD1)
  
  if(length(chr) > 1) stop("Please specify only one chromosome. ")
  n <- nrow(Y)
  p <- ncol(Y)
  
  genoprob <- pull.genoprob(cross, chr=chr)
  if(attr(cross,"class")[1] == "bc"){
    ngeno <- 2
  }else{
    ngeno <- 3
  }
  n.marker <- ncol(genoprob)/ngeno

  out1 <- scanone.mvn(cross=cross, Y=Y, chr=chr, addcov=addcov, intcov=intcov)
  LOD1 <- max(out1$lod)
  LOD1.pos <- which.max(out1$lod)

  ## find best qtl for each trait.
  p1 <- ncol(cross$pheno)
  cross$pheno <- data.frame(cross$pheno, Y)
  names(cross$pheno)[p1+(1:p)] <- colnames(Y)
  out <- scanone(cross, pheno=p1+(1:p), method="hk", chr=chr,
                 addcov=cbind(addcov,intcov), intcov=intcov)
  maxPOS <- apply(out[,-(1:2)],2,which.max)

  ## p qtl model: use best qtl of each trait 
  E <- matrix(NA, n, p)
  X <- cbind(rep(1, n), addcov, intcov)
  E <- .Call(stats:::C_Cdqrls, X, Y, tol)$residuals
  Sigma <- crossprod(E)
  L0 <- determinant(Sigma)$modulus

  E <- matrix(NA, n, p)
  for(i in 1:p){
    if(ngeno == 3){
      prob <- genoprob[,3*maxPOS[i]-2:1]
      X <- cbind(rep(1,n), addcov, intcov, prob, intcov*prob[,1], intcov*prob[,2])
    }else{
      prob <- genoprob[,2*maxPOS[i]-1]
      X <- cbind(rep(1,n), addcov, intcov, prob, intcov*prob)
    }
    E[,i] <- .Call(stats:::C_Cdqrls, X, Y[,i,drop=FALSE], tol)$residuals
  }
  Sigma <- crossprod(E)
  L1 <- determinant(Sigma)$modulus
  LODp <- n/2 * log10(exp(1)) * (L0 - L1)
  
  attributes(LODp) <-  NULL
  LODdiff <- LODp - LOD1
  names(maxPOS) <- NULL

  result <- list(LOD1 = LOD1, LODp = LODp, LODdiff = LODdiff,
                 LOD1.pos = LOD1.pos, maxPOS = maxPOS)
  return(result)
}



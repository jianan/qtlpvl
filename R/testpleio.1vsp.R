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
##' @param addcovar Additive covariates.
##' @param intcovar Interactive covariates.
##' @param n.simu number of simulations for p-value.
##' @param tol Tolerance value for the \code{qr} decomposition in
##' \code{lm} fitting.
##' @return a list of LOD1, LODp, LODdiff, LOD1.pos and maxPOS ... P-value
##' 
##' @export
##' @examples
##' data(hyper)
##' n <- 250
##' p <- 5
##' Y <- matrix(rnorm(n*p),n,p)
##' testpleio.1vsp(cross=hyper, Y=Y, chr="2")
testpleio.1vsp <- function(cross, Y, chr="6", addcovar=NULL, intcovar=NULL, n.simu=NA, tol=1e-7){

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

  out1 <- scanone.mvn(cross=cross, Y=Y, chr=chr, addcovar=addcovar, intcovar=intcovar)
  LOD1 <- max(out1$lod)
  LOD1.pos <- which.max(out1$lod)
  
  ## find best qtl for each trait.
  p1 <- ncol(cross$pheno)
  cross$pheno <- data.frame(cross$pheno, Y)
  if(!is.null(colnames(Y))) names(cross$pheno)[p1+(1:p)] <- colnames(Y)
  out <- scanone(cross, pheno.col=p1+(1:p), method="hk", chr=chr,
                 addcovar=cbind(addcovar,intcovar), intcovar=intcovar)
  maxPOS <- apply(out[,-(1:2)],2,which.max)

  ## p qtl model: use best qtl of each trait 
  E <- matrix(NA, n, p)
  X <- cbind(rep(1, n), addcovar, intcovar)
  E <- .Call(stats:::C_Cdqrls, X, Y, tol)$residuals
  Sigma <- crossprod(E)
  L0 <- determinant(Sigma)$modulus

  E <- matrix(NA, n, p)
  for(i in 1:p){
    if(ngeno == 3){
      prob <- genoprob[,3*maxPOS[i]-2:1]
      X <- cbind(rep(1,n), addcovar, intcovar, prob, intcovar*prob[,1], intcovar*prob[,2])
    }else{
      prob <- genoprob[,2*maxPOS[i]-1]
      X <- cbind(rep(1,n), addcovar, intcovar, prob, intcovar*prob)
    }
    E[,i] <- .Call(stats:::C_Cdqrls, X, Y[,i,drop=FALSE], tol)$residuals
  }
  Sigma <- crossprod(E)
  L1 <- determinant(Sigma)$modulus
  LODp <- n/2 * log10(exp(1)) * (L0 - L1)
  
  attributes(LODp) <-  NULL
  LODdiff <- LODp - LOD1
  names(maxPOS) <- NULL


  if(is.na(n.simu)){
    result <- list(LOD1 = LOD1, LODp = LODp, LODdiff = LODdiff,
                   LOD1.pos = LOD1.pos, maxPOS = maxPOS)
    return(result)
  } else{    ## simulation: parametric bootstrap.
    if(n.simu < 0) stop("n.simu should be a positive integer.")
    if(ngeno == 3){
      prob <- genoprob[,3*LOD1.pos-2:1]
      X <- cbind(rep(1,n), addcovar, intcovar, prob, intcovar*prob[,1], intcovar*prob[,2])
    }else{
      prob <- genoprob[,2*LOD1.pos-1]
      X <- cbind(rep(1,n), addcovar, intcovar, prob, intcovar*prob)
    }
    E.marker <- .Call(stats:::C_Cdqrls, X, Y, tol)$residuals

    Y.fit <- Y-E.marker
    Sigma <- cov(E.marker)
    Sigma.half <-  chol(Sigma)
    
    LODdiff.simu <- numeric(n.simu)

    for(i.simu in 1:n.simu){
      mat <- matrix(rnorm(p*n),p,n)      ## p*n
      mat <- crossprod(mat,Sigma.half)   ## n*p
      Y.simu <- Y.fit + mat
      result.i <- testpleio.1vsp(cross, Y.simu, chr=chr, addcovar=addcovar, intcovar=intcovar, n.simu=NA)
      LODdiff.simu[i.simu] <- result.i$LODdiff
    }
    pvalue <- mean(LODdiff > LODdiff.simu)
    result <- list(LOD1 = LOD1, LODp = LODp, LODdiff = LODdiff,
                   LOD1.pos = LOD1.pos, maxPOS = maxPOS,
                   pvalue=pvalue, LODdiff.simu = LODdiff.simu)
    return(result)
  }

}

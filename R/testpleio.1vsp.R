##' Test of pleiotrophic QTL vs close linked QTLs.
##'
##' Test if multitraits are controlled by one pleiotrophic QTL or
##' multiple close linked QTLs, allowing each trait to have its own
##' QTL. p-value is calculated from parametric bootstrap, which takes
##' the estimated parameters under null hypothesis and generate data
##' from them to get the null empirical distribution of test
##' statistic.
##'
##' @inheritParams testpleio.1vs2
##' @return a list of LOD1, LODp, LODdiff, LOD1.pos and maxPOS ... P-value
##' \item{LOD1}{LOD score of one dimensional joint mapping.}
##' \item{LODp}{LOD score for estimate of p QTL model, each trait
##'             is influenced by a separate QTL.}
##' \item{LODdiff}{Difference of best one QTL model and best p QTL model,
##'                LODdiff = max(LODp) - max(LOD1)}
##' \item{pvalue}{P-value from parametric bootstrap simulations.}
##'
##' @examples
##' data(fake.phenos)
##' data(listeria)
##' listeria <- calc.genoprob(listeria, step=1)
##' obj <- testpleio.1vsp(cross=listeria, Y=fake.phenos, chr=1, n.simu=100)
##' summary(obj)
##' plot(obj)
##'
##' @export
testpleio.1vsp <- function(cross, Y, chr="6", addcovar=NULL, intcovar=NULL, n.simu=1000, tol=1e-7){

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
  LOD1 <- out1$lod
  LOD1max <- max(out1$lod)

  ## find best qtl for each trait.
  p1 <- ncol(cross$pheno)
  cross$pheno <- data.frame(cross$pheno, Y)
  if(!is.null(colnames(Y))) names(cross$pheno)[p1+(1:p)] <- colnames(Y)
  out <- scanone(cross, pheno.col=p1+(1:p), method="hk", chr=chr,
                 addcovar=cbind(addcovar,intcovar), intcovar=intcovar)
  maxPOSind <- apply(out[,-(1:2)],2,which.max)
  maxPOS <- out$pos[apply(out[,-(1:2)],2,which.max)]
  maxLOD <- apply(out[,-(1:2)],2,max)

  ## p qtl model: use best qtl of each trait
  E <- matrix(NA, n, p)
  X <- cbind(rep(1, n), addcovar, intcovar)
  E <- lm.resid(X, Y)
  Sigma <- crossprod(E)
  L0 <- determinant(Sigma)$modulus

  E <- matrix(NA, n, p)
  for(i in 1:p){
    if(ngeno == 3){
      prob <- genoprob[,3*maxPOSind[i]-2:1]
      X <- cbind(rep(1,n), addcovar, intcovar, prob, intcovar*prob[,1], intcovar*prob[,2])
    }else{
      prob <- genoprob[,2*maxPOSind[i]-1]
      X <- cbind(rep(1,n), addcovar, intcovar, prob, intcovar*prob)
    }
    E[,i] <- lm.resid(X, Y[, i, drop=FALSE])
  }
  Sigma <- crossprod(E)
  L1 <- determinant(Sigma)$modulus
  LODp <- n/2 * log10(exp(1)) * (L0 - L1)

  attributes(LODp) <-  NULL
  LODdiff <- LODp - LOD1max

  map <- out1$pos
  map.marker <- unlist(pull.map(cross, chr))

  if(is.na(n.simu)){
    result <- list(LOD1=LOD1, LODp=LODp, LODdiff=LODdiff,
                   chr=chr, map=map, map.marker=map.marker,
                   maxLOD=maxLOD, maxPOS=maxPOS)
  } else{    ## simulation: parametric bootstrap.
    if(n.simu < 0) stop("n.simu should be a positive integer.")
    ind <- which.max(LOD1) ## index for LOD1 peak
    if(ngeno == 3){
      prob <- genoprob[,3*ind-2:1]
      X <- cbind(rep(1,n), addcovar, intcovar, prob, intcovar*prob[,1], intcovar*prob[,2])
    }else{
      prob <- genoprob[,2*ind-1]
      X <- cbind(rep(1,n), addcovar, intcovar, prob, intcovar*prob)
    }
    E.marker <- lm.resid(X, Y)

    Y.fit <- Y - E.marker
    Sigma <- cov(E.marker)
    Sigma.half <- chol(Sigma)

    LODdiff.simu <- numeric(n.simu)
    for(i.simu in 1:n.simu){
      mat <- matrix(rnorm(p*n),p,n)      ## p*n
      mat <- crossprod(mat,Sigma.half)   ## n*p
      Y.simu <- Y.fit + mat
      result.i <- testpleio.1vsp(cross, Y.simu, chr=chr,
                                 addcovar=addcovar, intcovar=intcovar, n.simu=NA)
      LODdiff.simu[i.simu] <- result.i$LODdiff
    }
    pvalue <- mean(LODdiff.simu > LODdiff - tol)
    result <- list(LOD1=LOD1, LODp=LODp, LODdiff=LODdiff,
                   chr=chr, map=map, map.marker=map.marker,
                   maxLOD=maxLOD, maxPOS=maxPOS, n.simu=n.simu,
                   pvalue=pvalue, LODdiff.simu=LODdiff.simu)
  }
  class(result) <- c("testpleio.1vsp", "list")
  return(result)
}

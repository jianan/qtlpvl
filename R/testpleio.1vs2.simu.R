testpleio.1vs2.simu <- function(cross, Y, E.marker, chr="6", addcov=NULL, intcov=NULL,
                                    int.method="bayes", search="fast",
                                    n.simu=NA, marker.l=NA, marker.r=NA){

  require(qtl) || stop("the required package 'qtl' is not installed. ")
  if(is.na(n.simu)) stop("n.simu should be specified")
  n <- nrow(Y);  p <- ncol(Y)
  
  Y.fit <- Y-E.marker
  Sigma <- cov(E.marker)
  Sigma.half <-  chol(Sigma)
  
  LODdiff.simu <- numeric(n.sim)
  L2inds <- matrix(NA,n.simu,2)
  Group <- matrix(NA,n.simu,p)

  for(i.simu in 1:n.simu){
    mat <- matrix(rnorm(p*n),p,n)      ## p*n
    mat <- crossprod(mat,Sigma.half)   ## n*p
    Y.simu <- Y.fit + mat
    result.i <- testpleio.1vs2(cross, Y.simu, chr=chr, addcov=addcov, intcov=intcov,
                                   int.method=int.method, search=search)
    L2inds[i.simu,] <- attr(result.i$LODdiff,"L2inds")
    LODdiff.simu[i.simu] <- result.i$LODdiff
    Group[i.simu,] <-result.i$Group
  }
  result <- list(LODdiff.simu=LODdiff.simu, Group=Group, L2inds=L2inds)  

  return(result)
}

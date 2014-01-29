##' Test of one pleiotrophic QTL vs two close linked QTLs.
##'
##' Test if multitriats are controlled by one pleiotrophic QTL or two
##' close linked QTLs. p-value is calculated from paramatric
##' bootstrap, which takes the estimated paramaters under null
##' hypothesis and generate data from them to get the null empirical
##' distribution of test statistic.
##'
##' @inheritParams scanone.mvn
##' @param region.l left bound
##' @param region.r right bound
##' @param int.method "bayes" or "1.5lod" method to calculated the
##' interval of interest if \code{region.l} and \code{region.r} is not
##' specified.
##' @param search searching method for two-QTL model.
##' @param RandomStart use random starting point for the two-QTL model
##' or not. default is \code{TRUE}.
##' @param n.simu number of simulations for p-value.
##' @param tol Tolerance value for the \code{qr} decomposition in
##' \code{lm} fitting.
##' @return a list of LODdiff, Group, P-value...
##' \item{LOD1}{One dimesional joint mapping.}
##' \item{LOD2}{LOD score for estimate of best two QTL model, each trait
##'             is influenced by the left or the right QTL.}
##' \item{LODdiff}{Difference of best one QTL model and best two QTL model,
##' LODdiff = max(LOD2) - max(LOD1)}
##' \item{pvalue}{P-value from parametric bootstrap simulations.}
##' @export
##' @examples
##' data(hyper)
##' n <- 250
##' p <- 5
##' Y <- matrix(rnorm(n*p),n,p)
##' hyper <- calc.genoprob(hyper)
##' summary(testpleio.1vs2(cross=hyper, Y=Y, chr="2", n.simu=100))

testpleio.1vs2 <- function(cross, Y, chr="6", addcovar=NULL, intcovar=NULL,
                           region.l=NA, region.r=NA, int.method="bayes",
                           search="fast", RandomStart=TRUE, n.simu=1000, tol=1e-7){

  if(length(chr) > 1) stop("Please specify only one chromosome. ")
  n <- nrow(Y); p <- ncol(Y)
  if(int.method!="bayes" && int.method!="1.5lod")
      stop("int.method need to be 'bayes' or '1.5lod'. ")
  
  p1 <- ncol(cross$pheno)
  cross$pheno <- data.frame(cross$pheno, Y)
  out <- scanone(cross, pheno.col=p1+(1:p), method="hk", chr=chr,
                 addcovar=cbind(addcovar,intcovar), intcovar=intcovar)

  if(is.na(region.l)){
    int <- matrix(NA,p,3)
    if(int.method=="1.5lod"){
      for(i in 1:p){
        int[i,] <- lodint(out[,c(1:2,2+i)],chr=chr)$pos
      }
    }else{  ##  if(int.method=="bayes")
      for(i in 1:p){
        int[i,] <- bayesint(out[,c(1:2,2+i)],chr=chr)$pos
      }
    }
    region.l <- min(int)
    region.r <- max(int)
  }
  
  map.marker <- unlist(pull.map(cross, chr))
  map.chr <- out$pos[out$chr==chr]
  if(sum(map.chr > region.l & map.chr < region.r)>0){
    rg <- range(which(map.chr > region.l & map.chr < region.r))
    marker.l <- rg[1]
    marker.r <- rg[2]
  }else{
    stop("confidence interval is too small.")
  }
  map.chr <- map.chr[marker.l:marker.r]
  
  ## order the traits by their mapped positions.
  pos <- out$pos
  maxPOS <- apply(out[,-(1:2)],2,which.max)
  o <- order(maxPOS)
  Y.orig <- Y
  Y <- Y[,o]
  maxPOS <- maxPOS[o]

  ## ---- scan.mvn ----
  genoprob <- pull.genoprob(cross, chr=chr)
  if(attr(cross,"class")[1] == "bc"){
    ngeno <- 2
    genoprob <- pull.genoprob(cross, chr=chr)[,(2*marker.l-1):(2*marker.r)]
  }else{
    ngeno <- 3
    genoprob <- pull.genoprob(cross, chr=chr)[,(3*marker.l-2):(3*marker.r)]
  }
  n.marker <- ncol(genoprob)/ngeno

  
  X <- cbind(rep(1,n), addcovar, intcovar)
  fit <- .Call(stats:::C_Cdqrls, X, Y, tol)
  Sigma <- crossprod(fit$residuals)
  L0 <- determinant(Sigma)$modulus  ## return log value
  
  L1 <- numeric(n.marker)
  E <- array(NA,dim=c(n,p,n.marker))
  Sigma.m <- array(NA,dim=c(p,p,n.marker)) ## save all the sigma matrix for later use
  for(i in 1:n.marker){
    if(ngeno == 3){
      prob <- genoprob[,(1:2)+3*(i-1)]
      X <- cbind(rep(1,n), prob[,1:2], addcovar, intcovar, intcovar*prob[,1], intcovar*prob[,2])
    }else{
      prob <- genoprob[,2*i-1]
      X <- cbind(rep(1,n), addcovar, intcovar, prob, intcovar*prob)
    }
    fit <- .Call(stats:::C_Cdqrls, X, Y, tol)
    E[,,i] <- fit$residuals
    Sigma.m[,,i] <- crossprod(E[,,i])
    L1[i] <- determinant(Sigma.m[,,i])$modulus  ## log value
  }

  LOD1 <- n/2*log10(exp(1))*(L0 - L1) ## scanone.mvn

  E.marker <- E[,,which.min(L1)]    ## residual matrix for the fitted model
  E.marker[,o] <- E.marker   ## change back to the original order...
  colnames(E.marker)[o] <- colnames(Y)
  rownames(E.marker) <- rownames(Y)

  if(length(unique(maxPOS))==1){    ## all QTLs mapped into same position...
    LODdiff <- 0
    attr(LODdiff,"LOD1lod") <- max(LOD1)  ## lod for the common QTL
    attr(LODdiff,"LOD1pos") <- map.chr[which.max(LOD1)]  ## pos for the common QTL
    attr(LODdiff,"LOD2lod") <- attr(LODdiff,"LOD1lod") ## pos for QTL1 and QTL2
    attr(LODdiff,"LOD2pos") <- rep(attr(LODdiff,"LOD1pos"),2) ## lod for QTL1 and QTL2
    Group <- rep(2,p)
    result <- list(E.marker=E.marker, LODdiff=LODdiff, Group=Group,
                   chr=chr, map.marker=map.marker, LOD1=LOD1)    
    class(result) <- c("testpleio.1vs2", "list")
    return(result)
  }


  ## ---- scan.2.mvn ----
  Sigma <- matrix(NA,p,p)
  L2mins <- rep(Inf,p-1) ## saves min of L2 for each cutting point.
  L2inds.trace <- matrix(NA,p-1,2) ## saves the index of minimum points of L2 for each cutting point.
  Group.trace <- matrix(2,p-1,p) ## saves the groupping vectore for each cutting point.
  
  ## searching:
  ## method 1. search all possible positions in the region.
  ##          [maybe only allow QTL1 on the left of the cut point and
  ##           QTL2 on the right]
  ## method 2. fix QTL2 and search for optimal QTL1, then fix QTL2 and search for optimal QTL2

  if(search=="complete"){    ## ---- do Complete Search in 2-dim ----
    for(i.cut in 1:(p-1)){
      if(maxPOS[i.cut] == maxPOS[i.cut+1]){
        next  ## skip when there are a group of traits mapped to the same pos
      }
      L2 <- matrix(NA,nrow=n.marker,ncol=n.marker)
      for(i in 1:n.marker){
        Sigma[1:i.cut, 1:i.cut] <- Sigma.m[1:i.cut, 1:i.cut, i]
        E1 <- E[,1:i.cut,i]
        for(j in 1:n.marker){
          Sigma[(i.cut+1):p,(i.cut+1):p] <- Sigma.m[(i.cut+1):p,(i.cut+1):p,j]
          E2 <- E[,(i.cut+1):p,j]
          tmp <- crossprod(E1,E2)
          Sigma[1:i.cut,(i.cut+1):p] <- tmp
          Sigma[(i.cut+1):p,1:i.cut] <- t(tmp)
          L2[i,j] <- determinant(Sigma)$modulus ## return log value
        }
      }
      L2mins[i.cut] <- min(L2,na.rm=TRUE)
      L2inds.trace[i.cut,] <- arrayInd(which.min(L2), .dim=dim(L2)) + marker.l - 1
      Group.trace[i.cut,o[1:i.cut]] <- 1
      if(i.cut==1 | L2mins[i.cut] == min(L2mins)){ ## saves the best cutting point
        L2.save <- L2 
      }
    }
    L2inds <- arrayInd(which.min(L2.save), .dim=dim(L2)) + marker.l - 1
  }
  
  if(search=="fast"){    ## Search 1-dim and then the other dim
    m.L1 <- which.min(L1)
    for(i.cut in 1:(p-1)){
      if(maxPOS[i.cut] == maxPOS[i.cut+1]){
        next  ## skip when there are a group of traits mapped to the same pos
      }
      L2 <- matrix(NA,nrow=n.marker,ncol=n.marker)
      i <- m.L1; j <- m.L1       ## start from peak of LOD1
      if(RandomStart){       ## start from a random position
        i <- sample(n.marker,1)
        j <- sample(n.marker,1)
      }
      while(TRUE){
        i.old <- i;  j.old <- j
        
        for(i in 1:n.marker){         ## fix j and search for i.min
          Sigma[1:i.cut, 1:i.cut] <- Sigma.m[1:i.cut, 1:i.cut, i]
          E1 <- E[,1:i.cut,i]
          Sigma[(i.cut+1):p,(i.cut+1):p] <- Sigma.m[(i.cut+1):p,(i.cut+1):p,j]
          E2 <- E[,(i.cut+1):p,j]
          tmp <- crossprod(E1,E2)
          Sigma[1:i.cut,(i.cut+1):p] <- tmp
          Sigma[(i.cut+1):p,1:i.cut] <- t(tmp)
          L2[i,j] <- determinant(Sigma)$modulus ## return log value
        }
        i <- which.min(L2[,j])
        
        for(j in 1:n.marker){        ## fix i and search for j.min
          Sigma[1:i.cut, 1:i.cut] <- Sigma.m[1:i.cut, 1:i.cut, i]
          E1 <- E[,1:i.cut,i]
          Sigma[(i.cut+1):p,(i.cut+1):p] <- Sigma.m[(i.cut+1):p,(i.cut+1):p,j]
          E2 <- E[,(i.cut+1):p,j]
          tmp <- crossprod(E1,E2)
          Sigma[1:i.cut,(i.cut+1):p] <- tmp
          Sigma[(i.cut+1):p,1:i.cut] <- t(tmp)
          L2[i,j] <- determinant(Sigma)$modulus ## return log value
        }
        j <- which.min(L2[i,])
        
        if(j == j.old & i==i.old){break}
      }
      L2mins[i.cut] <- min(L2,na.rm=TRUE)
      L2inds.trace[i.cut,] <- arrayInd(which.min(L2), .dim=dim(L2)) + marker.l - 1
      Group.trace[i.cut,o[1:i.cut]] <- 1
      if(i.cut == 1 | L2mins[i.cut] == min(L2mins)){  ## saves the best cutting point
        L2.save <- L2
      }
    }
    L2inds <- L2inds.trace[which.min(L2mins),]
  }
  
  LOD2 <- n/2*log10(exp(1))*(L0 - L2.save)
  
  LODdiff.trace <- -n/2*log10(exp(1))*(L2mins - min(L1)) ## LOD2-LOD1 for each cutting point.
  LODdiff <- max(LODdiff.trace)
  attr(LODdiff,"LOD1lod") <- max(LOD1)  ## lod for the common QTL
  attr(LODdiff,"LOD1pos") <- map.chr[which.max(LOD1)]  ## pos for the common QTL
  attr(LODdiff,"LOD2lod") <- max(LOD2,na.rm=TRUE) ## lod for QTL1 and QTL2
  attr(LODdiff,"LOD2pos") <- map.chr[L2inds] ## pos for QTL1 and QTL2
  ## attr(LODdiff,"L1inds") <- which.min(L1)+marker.l-1  ## index for the common QTL
  ## attr(LODdiff,"L2inds") <- L2inds ## index for QTL1 and QTL2
  ## attr(LODdiff,"region") <- map.chr[c(marker.l, marker.r)]

  Group <- rep(2,p)   ## a vectior of length p, with 1/-1 indicating grouped into right or left.'
  Group[o[1:which.min(L2mins)]] <- 1

  maxPOS <- out$pos[apply(out[,-(1:2)],2,which.max)]
  maxLOD <- apply(out[,-(1:2)],2,max)

  if(is.na(n.simu)){
    result <- list(LODdiff=LODdiff, Group=Group, chr=chr, map=map.chr,
                   maxPOS=maxPOS, maxLOD=maxLOD, LOD1=LOD1, LOD2=LOD2,
                   map.marker=map.marker,
                   ## Group.trace=Group.trace,
                   ## L2inds.trace=L2inds.trace, 
                   LODdiff.trace=LODdiff.trace)
  }else{    ## simulation: parametric bootstrap.
    if(n.simu < 0) stop("n.simu should be a positive integer.")
    
    Y.fit <- Y.orig - E.marker
    Sigma <- cov(E.marker)
    Sigma.half <- chol(Sigma)
    
    LODdiff.simu <- numeric(n.simu)
    LOD2pos <- matrix(NA,n.simu,2)
    Group.simu <- matrix(NA,n.simu,p)

    for(i.simu in 1:n.simu){
      mat <- matrix(rnorm(p*n),p,n)      ## p*n
      mat <- crossprod(mat,Sigma.half)   ## n*p
      Y.simu <- Y.fit + mat
      result.i <- testpleio.1vs2(cross, Y.simu, chr=chr, addcovar=addcovar, intcovar=intcovar,
                                 int.method=int.method, search=search, n.simu=NA)
      LOD2pos[i.simu,] <- attr(result.i$LODdiff,"LOD2pos")
      LODdiff.simu[i.simu] <- result.i$LODdiff
      Group.simu[i.simu,] <-result.i$Group
    }
    pvalue <- mean(LODdiff > LODdiff.simu)
    result <- list(LODdiff=LODdiff, Group=Group, chr=chr, map=map.chr,
                   maxPOS=maxPOS, maxLOD=maxLOD, LOD1=LOD1, LOD2=LOD2,
                   LODdiff.trace=LODdiff.trace,
                   map.marker=map.marker,
                   ## Group.trace=Group.trace,
                   ## L2inds.trace=L2inds.trace,
                   ## LOODdiff.simu=LODdiff.simu, Group.simu=Group.simu,
                   ## LOD2pos.simu=LOD2pos.simu,
                   pvalue=pvalue)
  }
  class(result) <- c("testpleio.1vs2", "list")
  return(result)
}


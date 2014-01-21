##' Test of one pleiotrophic QTL vs two close linked QTLs.
##'
##' Test if multitriats are controlled by one pleiotrophic QTL or two
##' close linked QTLs. p-value is calculated from paramatric
##' bootstrap, which takes the estimated paramaters under null
##' hypothesis and generate data from them to get the null empirical
##' distribution of test statistic.
##' 
##' @param Y matrix of multiple traits with rows equals to n, the
##' number of individuals.
##' @param cross An object of class \code{cross}. See
##' \code{read.cross} for details.
##' @param chr \code{chr, region.l, region.r} are used to specify the
##' interval of interest. see also \code{int.method}.
##' @param addcov Additive covariates.
##' @param intcov Interactive covariates.
##' @param region.l left bound
##' @param region.r right bound
##' @param int.method "bayes" or "1.5lod" method to calculated the
##' interval of interest if \code{region.l} and \code{region.r} is not
##' specified.
##' @param search searching method for two-QTL model.
##' @param RandomStart use random starting point for the two-QTL model
##' or not. default is \code{TRUE}.
##' @return a list of sth
##' @export
##' @examples
##' library(qtl)
##' data(hyper)
##' n <- 250
##' p <- 5
##' Y <- matrix(rnorm(n*p),n,p)
##' testpleio.1vs2(cross=hyper, Y=Y, chr="2")

testpleio.1vs2 <- function(cross, Y, chr="6", addcov=NULL, intcov=NULL,
                           region.l=NA, region.r=NA, int.method="bayes",
                           search="fast", RandomStart=TRUE){

  ## 1. scanone for each trait and order them by QTL.pos
  ## 2. LOD1 = LOD.scanone   
  ## 3. LOD2 = LOD for rightQTL and leftQTL
  ## 4. LODdiff = max(LOD2) - max(LOD1)
  
  require(qtl) || stop("the required package 'qtl' is not installed. ") 
  if(length(chr) > 1) stop("Please specify only one chromosome. ")
  n <- nrow(Y); p <- ncol(Y)
  if(int.method!="bayes" && int.method!="1.5lod") stop("int.method need to be 'bayes' or '1.5lod'. ")
  
  p1 <- ncol(cross$pheno)
  cross$pheno <- data.frame(cross$pheno, Y)
  out <- scanone(cross, pheno=p1+(1:p), method="hk", chr=chr,
                 addcov=cbind(addcov,intcov), intcov=intcov)

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
  
  map.chr <- out$pos[out$chr==chr]
  if(sum(map.chr > region.l & map.chr < region.r)>0){
    rg <- range(which(map.chr > region.l & map.chr < region.r))
    marker.l <- rg[1]
    marker.r <- rg[2]
  }else{
    stop("confidence interval is too small.")
  }
  
  ## order the traits by their mapped positions.
  maxPOS <- apply(out[,-(1:2)],2,which.max)
  o <- order(maxPOS)
  Y <- Y[,o]
  maxPOS <- maxPOS[o]

  ## ---- scan.mvn ----
  genoprob <- pull.genoprob(cross, chr=chr)[,(3*marker.l-2):(3*marker.r)]
  n.marker <- ncol(genoprob)/3
  
  L1 <- numeric(n.marker)
  E <- array(NA,dim=c(n,p,n.marker))
  Sigma.m <- array(NA,dim=c(p,p,n.marker)) ## save all the sigma matrix for later use
  for(i in 1:n.marker){
    prob <- genoprob[,(1:2)+3*(i-1)]
    X <- cbind(rep(1,n), prob[,1:2], addcov, intcov,
               intcov*prob[,1],intcov*prob[,2])
    fit <- .Call(stats:::C_Cdqrls, X, Y, 1e-7)
    E[,,i] <- fit$residuals
    Sigma.m[,,i] <- crossprod(E[,,i])
    L1[i] <- determinant(Sigma.m[,,i])$modulus  ## log value
  }

  ## residual matrix for the fitted model
  E.marker <- E[,,which.min(L1)] 
  ## reminber to change back to the original order...
  E.marker[,o] <- E.marker
  colnames(E.marker)[o] <- colnames(Y)
  rownames(E.marker) <- rownames(Y)


  if(length(unique(maxPOS))==1){    ## all QTLs mapped into same position...
    LODdiff <- 0
    attr(LODdiff,"L2inds") <- c(maxPOS[1],maxPOS[1])
    Group <- rep(2,p)
    result <- list(E.marker=E.marker, LODdiff=LODdiff, Group=Group)    
    return(result)
  }


  ## ---- scan.2.mvn ----
  Sigma <- matrix(NA,p,p)
  L2mins <- rep(Inf,p-1)
  L2inds.trace <- matrix(NA,p-1,2)
  Group.trace <- matrix(2,p-1,p)
  
  ## searching:
  ## method 1. search all possible positions in the region.
  ##          [maybe only allow QTL1 on the left of the cut point and
  ##           QTL2 on the right]
  ## method 2. fix QTL2 and search for optimal QTL1, then fix QTL2 and search for optimal QTL2

  if(search=="complete"){
    ## ---- do Complete Search in 2-dim ----
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
      ## arrayInd(which.min(L2), .dim=dim(L2))
      L2mins[i.cut] <- min(L2,na.rm=TRUE)
      L2inds.trace[i.cut,] <- arrayInd(which.min(L2), .dim=dim(L2))+ marker.l-1
      Group.trace[i.cut,o[1:i.cut]] <- 1
      if(i.cut==1 | L2mins[i.cut] == min(L2mins)){
        L2.save <- L2
      }
    }
    L2inds <- arrayInd(which.min(L2.save), .dim=dim(L2)) + marker.l -1
  }
  
  if(search=="fast"){
    ## Search 1-dim and then the other dim
    m.L1 <- which.min(L1)
    for(i.cut in 1:(p-1)){
      if(maxPOS[i.cut] == maxPOS[i.cut+1]){
        next  ## skip when there are a group of traits mapped to the same pos
      }
      L2 <- matrix(NA,nrow=n.marker,ncol=n.marker)
      ## start from peak of LOD1
      i <- m.L1; j <- m.L1
      ## start from a random position
      if(RandomStart){
        i <- sample(n.marker,1)
        j <- sample(n.marker,1)
      }
      while(TRUE){
        i.old <- i;  j.old <- j
        ## fix j and search for i.min
        for(i in 1:n.marker){
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
        
        ## fix i and search for j.min
        for(j in 1:n.marker){
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
      L2inds.trace[i.cut,] <- arrayInd(which.min(L2), .dim=dim(L2))+ marker.l-1
      Group.trace[i.cut,o[1:i.cut]] <- 1
      if(i.cut == 1 | L2mins[i.cut] == min(L2mins)){
        L2.save <- L2
      }
    }
    L2inds <- L2inds.trace[which.min(L2mins),]
  }
  
  
  LOD2maxs <- -n/2*log10(exp(1))*(L2mins - min(L1))
  LODdiff <- max(LOD2maxs)
  
  ind <- which.min(L2mins)

  attr(LODdiff,"L1inds") <- which.min(L1)+marker.l-1  ## pos for the common QTL
  attr(LODdiff,"L2inds") <- L2inds ## pos for QTL1 and QTL2
  attr(LODdiff,"L1pos") <- map.chr[which.min(L1)]  ## pos for the common QTL
  attr(LODdiff,"L2pos") <- map.chr[L2inds] ## pos for QTL1 and QTL2

  ## a vectior of length p, with 1/-1 indicating grouped into right or left.'
  Group <- rep(2,p)
  Group[o[1:ind]] <- 1
  ## attr(LODdiff,"Group") <- Group
  attr(LODdiff,"region") <- c(marker.l, marker.r)
  ## result <- list(E.marker=E.marker, LODdiff=LODdiff, Group=Group)
  
  ## return LOD1 and LOD2?
  X <- cbind(rep(1,n), addcov, intcov)
  fit <- .Call(stats:::C_Cdqrls, X, Y, 1e-7)
  Sigma <- crossprod(fit$residuals)
  L0 <- determinant(Sigma)$modulus  ## return log value
  LOD1 <- n/2*log10(exp(1))*(L0 - L1)
  LOD2 <- n/2*log10(exp(1))*(L0 - L2.save)

  ## return maxPOS and maxLOD for plot
  maxPOS <- apply(out[,-(1:2)],2,which.max)
  maxLOD <- apply(out[,-(1:2)],2,max)
  ## plot(maxPOS,maxLOD)
  ## plot(LOD1)
  ## image(LOD2)

  map <- map.chr[map.chr > region.l & map.chr < region.r]
  result <- list(E.marker=E.marker, LODdiff=LODdiff, Group=Group,
                 maxPOS=maxPOS, maxLOD=maxLOD, LOD1=LOD1, LOD2=LOD2,
                 LOD2maxs=LOD2maxs, map=map,
                 Group.trace=Group.trace, L2inds.trace=L2inds.trace)

  return(result)
}


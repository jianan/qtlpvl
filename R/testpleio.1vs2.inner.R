testpleio.1vs2.inner <- function(Y, maxPOS, genoprob, ngeno, addcovar, intcovar, 
                                 search, RandomStart, RandomCut, tol,
                                 in.simu=TRUE){

  if(in.simu & !RandomCut & length(unique(maxPOS))==1){    ## all QTLs mapped into same position...
    LODdiff <- 0
    return(LODdiff)
  }
  
  n <- nrow(Y)
  p <- ncol(Y)
  if(RandomCut) maxPOS <- 1:p
  
  X <- cbind(rep(1,n), addcovar, intcovar)
  Sigma <- crossprod(lm.resid(X, Y))
  L0 <- determinant(Sigma)$modulus  ## return log value
  
  n.marker <- ncol(genoprob)/ngeno
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
    E[,,i] <- lm.resid(X, Y)
    Sigma.m[,,i] <- crossprod(E[,,i])
    L1[i] <- determinant(Sigma.m[,,i])$modulus  ## log value
  }

  if(!in.simu & !RandomCut & length(unique(maxPOS))==1){    ## all QTLs mapped into same position...
    LODdiff <- 0
    LODdiff.trace <- rep(0, p)
    E.marker <- E[,,which.min(L1)]    ## residual matrix for the fitted model
    colnames(E.marker) <- colnames(Y)
    rownames(E.marker) <- rownames(Y)
    LOD1 <- n/2*log10(exp(1))*(L0 - L1) ## scanone.mvn
    LOD2 <- matrix(NA, n.marker, n.marker)
    result <- list(E.marker=E.marker, 
                   LODdiff.trace=LODdiff.trace,
                   LOD1=LOD1, LOD2=LOD2)    
    return(result)
  }

  ## ---- scan.2.mvn ----
  Sigma <- matrix(NA,p,p)
  L2mins <- rep(Inf,p-1) ## saves min of L2 for each cutting point.
  ## if(!in.simu){
  ##   L2inds.trace <- matrix(NA,p-1,2) ## saves the index of minimum points of L2 for each cutting point.
  ##   Group.trace <- matrix(2,p-1,p) ## saves the groupping vectore for each cutting point.
  ## }
  
  if(search=="complete"){    ## ---- do Complete Search in 2-dim ----
    for(i.cut in 1:(p-1)){
      if(!RandomCut & maxPOS[i.cut] == maxPOS[i.cut+1]){
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
      ## if(!in.simu){ ## save trace
      ##   L2inds.trace[i.cut,] <- arrayInd(which.min(L2), .dim=dim(L2)) + marker.l - 1
      ##   Group.trace[i.cut,o[1:i.cut]] <- 1
      ## }
      if(i.cut==1 | L2mins[i.cut] == min(L2mins)){ ## saves the best cutting point
        L2.save <- L2 
      }
    }
    ## L2inds <- arrayInd(which.min(L2.save), .dim=dim(L2)) + marker.l - 1
  }
  
  if(search=="fast"){    ## Search 1-dim and then the other dim
    m.L1 <- which.min(L1)
    for(i.cut in 1:(p-1)){
      if(!RandomCut & maxPOS[i.cut] == maxPOS[i.cut+1]){
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
        
        if(j == j.old & i==i.old) break
      }
      L2mins[i.cut] <- min(L2,na.rm=TRUE)
      ## if(!in.simu){ ## save trace
      ##   L2inds.trace[i.cut,] <- arrayInd(which.min(L2), .dim=dim(L2)) + marker.l - 1
      ##   Group.trace[i.cut,o[1:i.cut]] <- 1
      ## }
      if(i.cut == 1 | L2mins[i.cut] == min(L2mins)){  ## saves the best cutting point
        L2.save <- L2
      }
    }
    ## L2inds <- L2inds.trace[which.min(L2mins),]
  }

  if(!in.simu){
    E.marker <- E[,,which.min(L1)]    ## residual matrix for the fitted model
    colnames(E.marker) <- colnames(Y)
    rownames(E.marker) <- rownames(Y)
    LOD1 <- n/2*log10(exp(1))*(L0 - L1) ## scanone.mvn
    LOD2 <- n/2*log10(exp(1))*(L0 - L2.save)
    LODdiff.trace <- -n/2*log10(exp(1))*(L2mins - min(L1)) ## LOD2-LOD1 for each cutting point.
    result <- list(E.marker=E.marker, 
                   LODdiff.trace=LODdiff.trace,
                   LOD1=LOD1, LOD2=LOD2)    
    return(result)
  } else{
    LODdiff <- -n/2*log10(exp(1))*(min(L2mins) - min(L1)) 
    return(LODdiff)
  }

}

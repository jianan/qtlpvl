testpleio.1vs2.engine <- function(Y, maxPOS, genoprob, ngeno, addcovar, intcovar, 
                                  method, search.method, RandomStart, RandomCut, tol,
                                  in.simu=TRUE){

  if(in.simu & !RandomCut & length(unique(maxPOS))==1){    ## all QTLs mapped into same position...
    LODdiff <- 0
    return(LODdiff)
  }
  
  n <- nrow(Y)
  p <- ncol(Y)
  if(RandomCut) maxPOS <- 1:p
  
  X <- cbind(rep(1,n), addcovar, intcovar)
  E <- lm.resid(X, Y)
  
  L0 <- det_AtA(E)
  S0 <- crossprod(E)
  S0inv <- solve(S0)
  
  n.marker <- ncol(genoprob)/ngeno
  L1 <- numeric(n.marker)
  E <- array(NA,dim=c(n,p,n.marker))
  for(i in 1:n.marker){
    if(ngeno == 3){
      q <- 2 * p
      prob <- genoprob[,(1:2)+3*(i-1)]
      X <- cbind(rep(1,n), prob[,1:2], addcovar, intcovar, intcovar*prob[,1], intcovar*prob[,2])
    }else{
      q <- 1 * p
      prob <- genoprob[,2*i-1]
      X <- cbind(rep(1,n), addcovar, intcovar, prob, intcovar*prob)
    }
    df.res <- n-q-1
    E[,,i] <- E0 <- lm.resid(X, Y)
    L1[i] <- calc.L(E0, S0inv, q, df.res, method)
  }
  
  q <- q * 2 ## for 2QTL model
  df.res <- n-q-1
  
  if(!in.simu & !RandomCut & length(unique(maxPOS))==1){    ## all QTLs mapped into same position...
    LODdiff <- 0
    LODdiff.trace <- rep(0, p)
    E.marker <- E[,,which.min(L1)]    ## residual matrix for the fitted model
    colnames(E.marker) <- colnames(Y)
    rownames(E.marker) <- rownames(Y)
    if(method == "ML"){
      LOD1 <- n/2*log10(exp(1))*(L0 - L1) ## scanone.mvn
    }else{
      LOD1 <- - L1
    }
    LOD2 <- matrix(NA, n.marker, n.marker)
    result <- list(E.marker=E.marker, 
                   LODdiff.trace=LODdiff.trace,
                   LOD1=LOD1, LOD2=LOD2)    
    return(result)
  }

  ## ---- scan.2.mvn ----
  L2mins <- rep(Inf,p-1) ## saves min of L2 for each cutting point.
  
  if(search.method=="complete"){    ## ---- do Complete Search in 2-dim ----
    for(i.cut in 1:(p-1)){
      if(!RandomCut & maxPOS[i.cut] == maxPOS[i.cut+1]){
        next  ## skip when there are a group of traits mapped to the same pos
      }
      L2 <- matrix(NA,nrow=n.marker,ncol=n.marker)
      for(i in 1:n.marker){
        E1 <- E[,1:i.cut,i]
        for(j in 1:n.marker){
          E2 <- E[,(i.cut+1):p,j]
          E0 <- cbind(E1, E2)
          L2[i, j] <- calc.L(E0, S0inv, q, df.res, method)
        }
      }
      L2mins[i.cut] <- min(L2,na.rm=TRUE)
      if(i.cut==1 | L2mins[i.cut] == min(L2mins)){ ## saves the best cutting point
        L2.save <- L2 
      }
    }
  }
  
  if(search.method=="fast"){    ## Search 1-dim and then the other dim
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
          E1 <- E[,1:i.cut,i]
          E2 <- E[,(i.cut+1):p,j]
          E0 <- cbind(E1, E2)
          L2[i, j] <- calc.L(E0, S0inv, q, df.res, method)
        }
        i <- which.min(L2[,j])
        
        for(j in 1:n.marker){        ## fix i and search for j.min
          E1 <- E[,1:i.cut,i]
          E2 <- E[,(i.cut+1):p,j]
          E0 <- cbind(E1, E2)
          L2[i, j] <- calc.L(E0, S0inv, q, df.res, method)
        }
        j <- which.min(L2[i,])
        
        if(j == j.old & i==i.old) break
      }
      L2mins[i.cut] <- min(L2,na.rm=TRUE)
      if(i.cut == 1 | L2mins[i.cut] == min(L2mins)){  ## saves the best cutting point
        L2.save <- L2
      }
    }
  }

  if(!in.simu){
    E.marker <- E[,,which.min(L1)]    ## residual matrix for the fitted model
    colnames(E.marker) <- colnames(Y)
    rownames(E.marker) <- rownames(Y)
    if(method == "ML"){
      LOD1 <- n/2*log10(exp(1))*(L0 - L1) ## scanone.mvn
      LOD2 <- n/2*log10(exp(1))*(L0 - L2.save)
      LODdiff.trace <- -n/2*log10(exp(1))*(L2mins - min(L1)) ## LOD2-LOD1 for each cutting point.
    }else{
      LOD1 <- - L1
      LOD2 <- - L2.save
      ## LODdiff.trace <- - L2mins
      LODdiff.trace <- - L2mins + min(L1)
    }
    result <- list(E.marker=E.marker, 
                   LODdiff.trace=LODdiff.trace,
                   LOD1=LOD1, LOD2=LOD2)    
    return(result)
  }else{ ## in simulation, only need LODdiff.
    if(method == "ML"){
      LODdiff <- -n/2*log10(exp(1))*(min(L2mins) - min(L1)) 
    }else{
      ## LODdiff <- - min(L2mins)
      LODdiff <- - min(L2mins) + min(L1)
    }
    return(LODdiff)
  }
}

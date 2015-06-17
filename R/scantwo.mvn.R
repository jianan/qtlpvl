##' Multivariate normal genome 2d scan.
##'
##' Genome scan with two QTL model, with allowance for possible
##' covariates, using multivariate normal model for the multiple
##' traits.
##'
##' @inheritParams scanone.mvn
##' @export
scantwo.mvn <- function(cross, Y, chr=NULL, addcovar=NULL, intcovar=NULL,
                        method=c("maxlikelihood", "pillaitrace"),
                        tol=1e-7){

  method <- match.arg(method)
  ## checking inputs...
  if(class(cross)[2] != "cross") stop("cross need to be of class cross")
  n <- nrow(cross$pheno)
  if(missing(Y))  stop("Y needs to be speicfied")
  if(!missing(Y) & !is.matrix(Y)) stop("Y need to be a matrix.")
  if(!missing(Y) & is.matrix(Y) & n != nrow(Y)) stop("number of obs. in cross and Y not same.")
  p <- ncol(Y)
  if(p < 1) stop("Y should be a matrix with more than one columns")
  checkcov(intcovar, addcovar, n)

  if (!("prob" %in% names(cross[[c("geno", 1)]]))){
    warning("First running calc.genoprob.")
    cross <- calc.genoprob(cross)
  }

  if(missing(chr)){
    chr <- names(cross$geno)
  }else{
    cross <- subset(cross, chr)
  }

  ## todo-need to deal with X-chr
  ## chrclass <- sapply(cross$geno, class)
  genoprob <- pull.genoprob(cross, chr=chr, omit.first.prob=TRUE)
  ngeno <- ifelse(attr(cross, "class")[1] == "bc", 2, 3)
  m <- ncol(genoprob)/(ngeno-1)

  X0 <- cbind(rep(1, n), addcovar, intcovar)
  if(method == "maxlikelihood"){
    L0 <- det_AtA(lm.resid(X0, Y))
  }else if(method == "pillaitrace"){
    Sigma <- crossprod(lm.resid(X0, Y))
    S0inv <- solve(Sigma)
  }

  L1 <- matrix(NA, m, m)
  for(i in 1:m){
    if(ngeno==3) X1 <- genoprob[, 2*(i-1)+1:2] else X1 <- genoprob[, i]
    X <- cbind(rep(1, n), addcovar, intcovar, X1, interact(intcovar, X1))
    if(method == "maxlikelihood"){
      L1[i, i] <- det_AtA(lm.resid(X, Y))
    }else if(method == "pillaitrace"){
      L1[i, i] <- pillai(crossprod(lm.resid(X, Y)), S0inv)
    }
  }
  for(i in 1:(m-1)){
    if(ngeno==3) X1 <- genoprob[, 2*(i-1)+1:2] else X1 <- genoprob[, i]
    for(j in (i+1):m){
      if(ngeno==3) X2 <- genoprob[, 2*(j-1)+1:2] else X2 <- genoprob[, j]
      Xa <- cbind(rep(1, n), addcovar, intcovar, X1, X2,
                  interact(intcovar, X1), interact(intcovar, X2))
      Xi <- interact(X1, X2)
      Xf <- cbind(Xa, Xi, interact(intcovar, Xi))
      if(method == "maxlikelihood"){
        L1[i, j] <- det_AtA(lm.resid(Xa, Y)) ## upper tri, add model
        L1[j, i] <- det_AtA(lm.resid(Xf, Y)) ## lower tri, full model
      }else if(method == "pillaitrace"){
        L1[i, j] <- pillai(crossprod(lm.resid(Xa, Y)), S0inv)
        L1[j, i] <- pillai(crossprod(lm.resid(Xf, Y)), S0inv)
      }
    }
  }

  if(method == "maxlikelihood"){
    LOD <- n/2 * log10(exp(1)) * (L0 - L1)
  }else if(method == "pillaitrace"){
    LOD <- (L1 - p/2)/p
  }

  ## out <- NULL
  ## for(i in chr){
  ##   map <- attr(cross[[c("geno", i, "prob")]], "map")
  ##   ind <- substr(names(map), 1, 3) == "loc"
  ##   names(map)[ind] <- paste("c", i, ".", names(map)[ind], sep="")
  ##   out <- rbind(out, data.frame(chr=i, pos=map))
  ## }
  ## out$lod <- LOD
  ## class(out) <- c("scanone", "data.frame")
  ## ## attr(out, "method") <- method
  ## ## attr(out, "type") <- type
  ## ## attr(out, "model") <- model
  ## return(out)
  return(LOD)
}

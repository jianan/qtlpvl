##' Multivariate normal genome scan with a single QTL model.
##' 
##' Genome scan with a single QTL model, with allowance for possible
##' covariates, using multivariate normal model for the multiple
##' traits.
##'
##' @param cross An object of class \code{cross}. See
##' \code{read.cross} for details.
##' @param Y Matrix of multiple traits with n rows.
##' @param chr Optional vector indicating the chromosomes for which
##' LOD scores should be calculated.  This should be a vector of
##' character strings referring to chromosomes by name; numeric values
##' are converted to strings.
##' @param addcovar Additive covariates.
##' @param intcovar Interactive covariates.
##' @param tol Tolerance value for the \code{qr} decomposition in
##' \code{lm} fitting.
##' @return A data.frame whose first column contains the chromosome
##' IDs, second column contains cM positions, third column contains
##' the LOD scores.
##' @keywords QTL
##' @seealso qtl::scanone
##'
##' @examples
##' set.seed(92950640)
##' data(listeria)
##' listeria <- calc.genoprob(listeria,step=1)
##' n <- nind(listeria)
##' chr <- "1"
##' geno <- pull.geno(listeria, chr=chr)
##' genotype1 <- geno[,7]
##' genotype2 <- geno[,10]
##' p <- 10
##' p1 <- floor(p/2)
##' G1 <- matrix(genotype1, n, p1)
##' G2 <- matrix(genotype2, n, p-p1)
##' G2[G2==3] <- 2
##' G <- cbind(G1, G2*(-2))
##' Y <- matrix(rnorm(n*p),n,p)
##' Y <- Y + G
##' out <- scanone.mvn(listeria, Y, chr)
##' summary(out)
##' plot(out)
##'
##' @export
scanone.mvn <- function(cross, Y, chr=NULL, addcovar=NULL, intcovar=NULL, tol=1e-7){

  ## checking inputs...
  if(class(cross)[2] != "cross") stop("cross need to be of class cross")
  n <- nrow(cross$pheno)
  if(missing(Y))  stop("Y needs to be speicfied")
  if(!missing(Y) & !is.matrix(Y)) stop("Y need to be a matrix.")
  if(!missing(Y) & is.matrix(Y) & n != nrow(Y)) stop("number of obs. in cross and Y not same.")
  p <- ncol(Y)
  if(p < 1) stop("Y should be a matrix with more than one columns")
  checkcov(intcovar,addcovar,n)
  
  if (!("prob" %in% names(cross[[c("geno",1)]]))){
    warning("First running calc.genoprob.")
    cross <- calc.genoprob(cross)
  }

  if(missing(chr)){
    chr <- names(cross$geno)
  }else{
    cross <- subset(cross, chr)
  }

  ## todo-need to deal with X-chr
  chrclass <- sapply(cross$geno,class)
  genoprob <- pull.genoprob(cross, chr=chr)
  if(attr(cross,"class")[1] == "bc"){
    ngeno <- 2
  }else{
    ngeno <- 3
  }
  m <- ncol(genoprob)/ngeno
  
  E <- matrix(NA, n, p)
  X <- cbind(rep(1, n), addcovar, intcovar)
  E <- lm.resid(X, Y)
  Sigma <- crossprod(E)
  L0 <- determinant(Sigma)$modulus
  
  L1 <- numeric(m)
  for(i in 1:m){
    if(ngeno == 3){
      prob <- genoprob[,3*i-2:1]
      X <- cbind(rep(1,n), addcovar, intcovar, prob, intcovar*prob[,1], intcovar*prob[,2])
    }else{
      prob <- genoprob[,2*i-1]
      X <- cbind(rep(1,n), addcovar, intcovar, prob, intcovar*prob)
    }
    E <- lm.resid(X, Y)
    Sigma <- crossprod(E)
    L1[i] <- determinant(Sigma)$modulus
  }
  LOD <- n/2 * log10(exp(1)) * (L0 - L1)

  out <- NULL
  for(i in chr){
    map <- attr(cross[[c("geno",i,"prob")]], "map")
    ind <- substr(names(map),1,3) == "loc"
    names(map)[ind] <- paste("c", i, ".", names(map)[ind], sep="")
    out <- rbind(out, data.frame(chr=i, pos=map))
  }
  out$lod <- LOD
  class(out) <- c("scanone", "data.frame")
  ## attr(out, "method") <- method
  ## attr(out, "type") <- type
  ## attr(out, "model") <- model
  return(out)
}

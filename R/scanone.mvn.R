##' Multivariate normal genome scan with a single QTL model.
##'
##' Genome scan with a single QTL model, with allowance for possible
##' covariates, using multivariate normal model for the multiple
##' traits.
##'
##' @param cross An object of class \code{cross}. See
##' \code{read.cross} for details.
##' @param Y Matrix of multiple traits, with samples in the row,
##' traits in the column.
##' @param chr Optional vector indicating the chromosomes for which
##' LOD scores should be calculated. This should be a vector of
##' character strings referring to chromosomes by name; numeric values
##' are converted to strings.
##' @param addcovar Additive covariates, samples in the row.
##' @param intcovar Interactive covariates (interact with QTL genotype).
##' @param method choose from "ML" for maximum likelihood method for
##' LOD score, "Pillai", "Wilks", "Hotelling-Lawley", or "Roy" for the
##' corresponding MANOVA methods.
##' @param tol Tolerance value for the \code{qr} decomposition in
##' \code{lm} fitting.
##' @return A data.frame whose first column contains the chromosome
##' IDs, second column contains cM positions, third column contains
##' the LOD scores.
##' @keywords QTL
##' @seealso qtl::scanone
##'
##' @examples
##' data(fake.phenos)
##' data(listeria)
##' listeria <- calc.genoprob(listeria, step=1)
##' out <- scanone.mvn(Y=fake.phenos, cross=listeria, chr=1)
##' summary(out)
##' plot(out)
##'
##' @export
scanone.mvn <- function(cross, Y, chr=NULL, addcovar=NULL, intcovar=NULL,
                        method = c("ML", "Pillai", "Wilks", "Hotelling-Lawley", "Roy"),
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
  checkcov(intcovar,addcovar,n)

  ## if there is missing value in Y, remove thoes individuals.
  if(any(is.na(Y))){
    id <- which(!apply(is.na(Y), 1, any))
    Y <- Y[id, ]
    cross <- subset(cross, ind=id)
    addcovar <- addcovar[id, ]
    intcovar <- intcovar[id, ]
    n <- nrow(Y)
  }

  if (!("prob" %in% names(cross[[c("geno",1)]]))){
    warning("First running calc.genoprob.")
    cross <- calc.genoprob(cross)
  }

  if(missing(chr)){
    chr <- names(cross$geno)
    chr <- chr[chr!="X"] ## not deal with X chr for now.
  }
  cross <- subset(cross, chr=chr)

  ## todo-need to deal with X-chr
  chrclass <- sapply(cross$geno,class)
  genoprob <- pull.genoprob(cross, chr=chr)
  if(attr(cross,"class")[1] == "bc"){
    ngeno <- 2
  }else{
    ngeno <- 3
  }
  m <- ncol(genoprob)/ngeno

  X <- cbind(rep(1, n), addcovar, intcovar)
  E <- lm.resid(X, Y)
  L0 <- det_AtA(E)
  S0 <- crossprod(E)
  S0inv <- solve(S0)

  L1 <- numeric(m)
  for(i in 1:m){
    if(ngeno == 3){
      q <- 2 * p
      prob <- genoprob[,3*i-2:1]
      X <- cbind(rep(1,n), addcovar, intcovar, prob, intcovar*prob[,1], intcovar*prob[,2])
    }else{
      q <- 1 * p
      prob <- genoprob[,2*i-1]
      X <- cbind(rep(1,n), addcovar, intcovar, prob, intcovar*prob)
    }
    df.res <- n-q-1
    E <- lm.resid(X, Y)
    L1[i] <- calc.L(E, S0inv, q, df.res, method)
  }

  if(method=="ML"){
    LOD <- n/2 * log10(exp(1)) * (L0 - L1)
  }else{
    LOD <- - L1
  }

  out <- NULL
  for(i in chr){
    map <- attr(cross[[c("geno",i,"prob")]], "map")
    ind <- substr(names(map),1,3) == "loc"
    names(map)[ind] <- paste("c", i, ".", names(map)[ind], sep="")
    out <- rbind(out, data.frame(chr=i, pos=map))
  }
  out$lod <- LOD
  class(out) <- c("scanone", "data.frame")
  attr(out, "method") <- method
  return(out)
}

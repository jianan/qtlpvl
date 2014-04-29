##' Calculate residual matrix.
##'
##' Run lm(Y~X) with specific method through Rcpp and return a
##' residual matrix.
##'
##' @param X A model matrix
##' @param Y The response matrix
##' @param method 'llt' for the LLT Cholesky, 'qr' for the
##' column-pivoted QR decomposition, 'svd' for the Jacobi singular
##' value decomposition (SVD)
##' @return The residual matrix.
##' @export
lm.resid <- function(X, Y, method=c("llt", "qr", "svd")){
  stopifnot(is.matrix(X))
  stopifnot(is.matrix(Y))
  stopifnot(nrow(X)==nrow(Y))
  method <- match.arg(method)
  if(method=="llt"){
    res <- lm_resid_llt(X, Y)
  }else if(method=="qr"){
    res <- lm_resid_qr(X, Y)
  }else if(method=="svd"){
    res <- lm_resid_svd(X, Y)
  }
  dimnames(res) <- dimnames(Y)
  return(res)
}

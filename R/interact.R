##' computes a matrix represents the interaction of two given matrices.
##' @param x Matrix
##' @param y Matrix
##' @return Matrix of interaction of x and y
##' @export
interact <- function(x,y){
  if(is.null(x)|is.null(y)) return()
  if(!is.matrix(x)) x <- as.matrix(x)
  if(!is.matrix(y)) y <- as.matrix(y)
  stopifnot(nrow(x)==nrow(y))
  px <- ncol(x)
  py <- ncol(y)
  res <- matrix(NA, nrow(x), px*py)
  for(i in 1:px)
      for(j in 1:px)
          res[, (i-1)*px+j] <- x[, i] * y[, j]
  return(res)
}

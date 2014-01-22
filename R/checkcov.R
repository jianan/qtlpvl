checkcov <- function(intcovar, addcovar, n){
  ## From karl's qtl package, util.R
  if(!is.null(addcovar)) {
    if(!is.matrix(addcovar)) {
      if(!is.numeric(as.matrix(addcovar)))
          stop("addcovar should be numeric")
      if(is.vector(addcovar) || is.data.frame(addcovar))
          addcovar <- as.matrix(addcovar)
      else stop("addcovar should be a matrix")
    }
    if(!all(apply(addcovar,2,is.numeric)))
        stop("All columns of addcovar must be numeric")
    if( nrow(addcovar) != n) {
      stop("Number of rows in additive covariates is incorrect")
    }
  }
  if(!is.null(intcovar)) {
    if(!is.matrix(intcovar)) {
      if(!is.numeric(as.matrix(intcovar)))
          stop("intcovar should be a numeric")
      if(is.vector(intcovar) || is.data.frame(intcovar))
          intcovar <- as.matrix(intcovar)
      else stop("intcovar should be a matrix")
    }
    if(!all(apply(intcovar,2,is.numeric)))
        stop("All columns of intcovar must be numeric")
    if(nrow(intcovar)[1] != n) {
      stop("The length of interacting covariates is incorrect!")
    }
  }
}

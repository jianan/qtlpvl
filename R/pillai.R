pillai <- function(S1, S0inv){
  A <- crossprod(S1, S0inv)
  lambda <- eigen(A)$values
  if(is.complex(lambda)){
    if(all(abs(Im(lambda)) < 1e-7)) lambda <- Re(lambda)
    else return(NA)
  } 
  return(sum(1/(1+lambda)))
}

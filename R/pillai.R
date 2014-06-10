pillai <- function(S1, S0inv){
  A <- crossprod(S1, S0inv)
  lambda <- eigen(A)$values
  return(sum(1/(1+lambda)))
}

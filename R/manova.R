## From R-base, manova.

Pillai <- function (eig, q, df.res)
{
  test <- sum(eig/(1 + eig))
  p <- length(eig)
  s <- min(p, q)
  n <- 0.5 * (df.res - p - 1)
  m <- 0.5 * (abs(p - q) - 1)
  tmp1 <- 2 * m + s + 1
  tmp2 <- 2 * n + s + 1
  c(test, (tmp2/tmp1 * test)/(s - test), s * tmp1, s * tmp2)
}

Wilks <- function (eig, q, df.res)
{
  test <- prod(1/(1 + eig))
  p <- length(eig)
  tmp1 <- df.res - 0.5 * (p - q + 1)
  tmp2 <- (p * q - 2)/4
  tmp3 <- p^2 + q^2 - 5
  tmp3 <- if (tmp3 > 0)
      sqrt(((p * q)^2 - 4)/tmp3)
  else 1
  c(test, ((test^(-1/tmp3) - 1) * (tmp1 * tmp3 - 2 * tmp2))/p/q,
    p * q, tmp1 * tmp3 - 2 * tmp2)
}

HL <- function (eig, q, df.res)
{
  test <- sum(eig)
  p <- length(eig)
  m <- 0.5 * (abs(p - q) - 1)
  n <- 0.5 * (df.res - p - 1)
  s <- min(p, q)
  tmp1 <- 2 * m + s + 1
  tmp2 <- 2 * (s * n + 1)
  c(test, (tmp2 * test)/s/s/tmp1, s * tmp1, tmp2)
}

Roy <- function (eig, q, df.res)
{
  p <- length(eig)
  test <- max(eig)
  tmp1 <- max(p, q)
  tmp2 <- df.res - tmp1 + q
  c(test, (tmp2 * test)/tmp1, tmp1, tmp2)
}


calc.L <- function(E, S0inv, q, df.res, method){
  ## clac L1 or L2 used in LOD scores, smaller is better.
  if(method != "ML"){
    S1 <- crossprod(E)
    A <- crossprod(S1, S0inv)
    d <- Re(eigen(A, symmetric = FALSE)$values)
    eig <- 1/d - 1
  }
  res <- switch(method,
                ML = det_AtA(E),
                Pillai = - Pillai(eig, q, df.res)[1],
                Wilks = log10( Wilks(eig, q, df.res)[1]),
                `Hotelling-Lawley` = - HL(eig, q, df.res)[1],
                Roy = - Roy(eig, q, df.res)[1])
  res
}

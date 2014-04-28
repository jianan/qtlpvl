context("tests of lm in Rcpp")

test_that("tests of lm_resid, lm_resid_cov, lm_resid_cov_det",{
  set.seed(912803290)
  n <- 500
  p <- 50

  Y <- matrix(rnorm(n*p), n, p)
  X <- matrix(rnorm(n*2), n, 2)
  
  E.r <- lm.fit(X, Y)$residuals
  E.c <- lm.resid(X, Y)
  expect_equal(E.r, E.c)

  EtE.r <- crossprod(E.r)
  EtE.c <- lm_resid_cov(X, Y)
  expect_equal(EtE.r, EtE.c)

  L2.r <- determinant(EtE.r)$modulus
  L2.c <- lm_resid_cov_det(X ,Y) 
  expect_equal(L2.r, L2.c,  check.attributes=FALSE)
})


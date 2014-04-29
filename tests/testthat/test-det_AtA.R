context("tests of det in Rcpp")

test_that("tests of det",{
  set.seed(74695464)
  n <- 500
  p <- 50

  X <- matrix(rnorm(n*2), n, 2)
  
  det.r <- determinant(crossprod(X))$modulus
  det.c <- det_AtA(X)
  expect_equal(det.r, det.c,  check.attributes=FALSE)
})


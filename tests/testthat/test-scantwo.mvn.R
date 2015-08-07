library(qtl)
data(hyper)
hyper <- calc.genoprob(hyper,step=2)

context("tests on inputs")

test_that("tests for variable Y",{
  set.seed(12345)
  cross <- hyper
  Y <- rnorm(250)
  expect_that(scantwo.mvn(cross=cross, Y=Y),
              throws_error("Y need to be a matrix."))

  n <- 150
  p <- 5
  Y <- matrix(rnorm(n*p),n,p)
  expect_that(scanone.mvn(cross=cross, Y=Y),
              throws_error("number of obs. in cross and Y not same."))
})


context("tests on outputs")

test_that("tests for LOD score",{
  set.seed(12345)
  cross <- hyper
  n <- 250
  p <- 5
  Y <- matrix(rnorm(n*p),n,p)
  chr <- 1

  LOD1 <- scanone.mvn(cross=cross, Y=Y, chr=chr)$lod
  LOD2 <- scantwo.mvn(cross=cross, Y=Y, chr=chr)
  expect_equal(LOD1, diag(LOD2))
})

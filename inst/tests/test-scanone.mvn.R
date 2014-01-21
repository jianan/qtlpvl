
library(qtl)
data(hyper)
hyper <- calc.genoprob(hyper,step=2)

context("tests on inputs")

test_that("tests for variable Y",{
  set.seed(12345)
  cross <- hyper
  Y <- rnorm(250)
  expect_that(scanone.mvn(cross=cross, Y=Y),
              throws_error("Y need to be a matrix."))

  n <- 150
  p <- 5
  Y <- matrix(rnorm(n*p),n,p)
  expect_that(scanone.mvn(cross=cross, Y=Y),
              throws_error("number of obs. in cross and Y not same."))
})


context("tests on outputs")

test_that("tests for chr ",{
  set.seed(12345)
  cross <- hyper
  n <- 250
  p <- 5
  Y <- matrix(rnorm(n*p),n,p)
  chr <- c(1:10,"X")

  out <- scanone.mvn(cross=cross, Y=Y, chr=chr)
  out.chr <- unique(out$chr)
  out.chr <- as.character(out.chr)
  expect_identical(out.chr, chr)
})



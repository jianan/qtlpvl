context("tests of classifications")

test_that("tests of lda",{
  set.seed(96308511)

  library(qtl)
  data(listeria)
  n <- nind(listeria)
  listeria <- calc.genoprob(listeria, step=5, stepwidth="max")
  Xgp <- pull.genoprob(listeria, chr=c(1, 2), omit.first.prob=TRUE)
  gp <- Xgp[, grep("D1M291", colnames(Xgp))]
  G <- gp[,1] * 10 + gp[,2] * 20
  p <- 10
  Y <- matrix(rnorm(n*p),n,p)
  Y <- Y + G
  colnames(Y) <- paste0("pheno", 1:p)

  group <- group.train.test(listeria, Y, chr=1, region.l=75, region.r=95)
  data.train <- group$data.train
  data.test <- group$data.test
  class.train <- group$geno.train
  genotype <- group$geno.test
  map <- group$map

  fit <- classification(data.train, data.test, class.train, method="LDA")
  pred.test <- fit$pred.test
  pred.score <- fit$pred.score
  sca <- fit$sca
  error.train <- fit$error.train

  layout(matrix(c(1,2,3,3),2,2))
  mosaic.plot(genotype, map, pred.test, pred.score, main="", weighted = FALSE, label=FALSE)
  propn.plot(genotype, map, pred.test, pred.score, main="", weighted = FALSE)
  plotlda(data.train, data.test, class.train, pred.test, error.train, sca, main="main")

  expect_equal(error.train, 0)
  expect_identical(pred.test ,
            structure(c(2L, 2L, 1L, 2L, 2L, 1L, 1L, 3L, 3L, 2L, 3L, 2L, 1L,
                        2L, 2L, 3L, 3L, 2L, 1L, 1L, 2L, 2L, 2L, 2L, 3L, 1L,
                        2L, 2L, 2L, 3L, 2L),
                      .Label = c("1", "2", "3"), class = "factor",
                      .Names = c("16", "18", "20", "24", "27", "31", "32", "35", "36", "42",
                          "44", "48", "50", "63", "64", "68", "86", "88", "90", "91", "92",
                          "93", "95", "97", "102", "103", "111", "112", "116", "117", "120")))

})

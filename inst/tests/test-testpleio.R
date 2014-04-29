context("tests of testpleio.1vs2 and testpleio.1vsp")

## generate data.
library(qtl)
set.seed(92950640)
data(listeria)
listeria <- calc.genoprob(listeria,step=1)
n <- nind(listeria)
chr <- "1"
geno <- pull.geno(listeria, chr=chr)
genotype1 <- geno[,7]
genotype2 <- geno[,10]
(pos <- unlist(pull.map(listeria, chr=chr))[c(7,10)])
p <- 10
p1 <- floor(p/2)
G1 <- matrix(genotype1, n, p1)
G2 <- matrix(genotype2, n, p-p1)
G2[G2==3] <- 2
G <- cbind(G1, G2*(-2))
Y <- matrix(rnorm(n*p),n,p)
Y <- Y + G

library(devtools)
load_all("~/Project/R-pkg/qtlpvl")

test_that("tests of testpleio",{

  n.simu <- 10
  region.l <- 60
  region.r <- 90
  RandomCut <- TRUE
  RandomStart <- TRUE
  int.method <- "bayes"
  search.method <- "complete"
  simu.method <- "permutation"
  obj <- testpleio.1vs2(listeria, Y, chr, n.simu=n.simu,
                        region.l=region.l, region.r=region.r, 
                        search.method=search.method,
                        RandomCut=RandomCut,
                        RandomStart=RandomStart,
                        simu.method=simu.method)
  
  expect_false(any(is.na(obj$LOD2)))
  expect_true(obj$pvalue >= 0)
  
  n.simu <- 2
  for(RandomCut in c(TRUE, FALSE))
      for(RandomStart in c(TRUE, FALSE))
          for(int.method in c("bayes", "1.5lod"))
              for(search.method in c("fast", "complete"))
                  for(simu.method in c("permutation", "parametric"))
                      {obj <- testpleio.1vs2(listeria, Y, chr, n.simu=n.simu,
                                             int.method=int.method,
                                             search.method=search.method,
                                             RandomCut=RandomCut,
                                             RandomStart=RandomStart,
                                             simu.method=simu.method)
                     }

  obj.1vsp <- testpleio.1vsp(listeria, Y, chr, n.simu=n.simu)
  expect_true(obj$LODdiff <= obj.1vsp$LODdiff)
  expect_equal(obj$maxLOD, obj.1vsp$maxLOD, check.attributes=FALSE)
  expect_equal(obj$maxPOS, obj.1vsp$maxPOS, check.attributes=FALSE)
  expect_true(obj.1vsp$pvalue >= 0)
})

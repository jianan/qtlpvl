
library("devtools")
library("roxygen2")
library("knitr")
library("testthat")

pkg <- "../qtlpvl"

document(pkg)
load_all(pkg)
test(pkg)
run_examples(pkg)

install(pkg)
## help(package="qtlpvl")
## example(scanone.mvn)
## example(testpleio.1vs2)
## example(testpleio.1vsp)
## example(plotGenetpattern)
## example(plotLODsign)



## ------------ examples --------------

data(hyper)
hyper <- calc.genoprob(hyper)
geno <- pull.geno(hyper, chr="1")
genotype1 <- geno[,6]
genotype2 <- geno[,12]
n <- length(genotype1)
p <- 10
p1 <- floor(p/2)
G1 <- matrix(genotype1, n, p1)
G2 <- -matrix(genotype2, n, p-p1)
G <- cbind(G1, G2)
Y <- matrix(rnorm(n*p, sd=0.5), n, p)
Y <- Y + G

obj <- testpleio.1vs2(cross=hyper, Y=Y, chr="1", n.simu=100)
summary(obj)
plot(obj)

obj <- testpleio.1vsp(cross=hyper, Y=Y, chr="1", n.simu=100)
summary(obj)
plot(obj)

#####
data(listeria)
geno <- pull.geno(listeria, chr="1")
genotype1 <- geno[,7]
genotype2 <- geno[,10]
n <- length(genotype1)
p <- 100
p1 <- floor(p/2)
G1 <- matrix(genotype1, n, p1)
G2 <- matrix(genotype2, n, p-p1)
G2[G2==3] <- 2
G2 <- -G2
G <- cbind(G1, G2)
Y <- matrix(rnorm(n*p,sd=0.5),n,p)
Y <- Y + G
plotGenetpattern(Y, genotype1)

####
data(listeria)
listeria <- calc.genoprob(listeria)
geno <- pull.geno(listeria, chr="1")
genotype1 <- geno[,7]
genotype2 <- geno[,10]
n <- length(genotype1)
p <- 100
p1 <- floor(p/2)
G1 <- matrix(genotype1, n, p1)
G2 <- -matrix(genotype2, n, p-p1)
G <- cbind(G1, G2)
Y <- matrix(rnorm(n*p,sd=0.5),n,p)
Y <- Y + G
plotLODsign(listeria, Y, chr="1")

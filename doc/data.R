library(qtl)
set.seed(92950640)
data(listeria)
listeria <- calc.genoprob(listeria, step=1)
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
phenoname <- paste0("phenos", 1:p)
colnames(Y) <- phenoname

map <- pull.map(listeria)
chrs <- sample(2:19, size=p,replace=TRUE)
cM <- runif(n=p, min=0, max=sapply(map, max)[chrs])
probepos <- data.frame(chr=chrs,cM=cM)
rownames(probepos) <- phenoname
o <- order(probepos$chr, probepos$cM)
probepos <- probepos[o, ]

save(listeria, Y, n, p, chr, pos, probepos, file="../data/fake.phenos.RData")

##' find the peak of scanone.mvn and the QTL genotype.
##'
##' For a list of hotspots, find the single QTL genotype and the mice
##' with no recombinant events in a 10cM region around the single QTL.
##'
##' @inheritParams add.effect
##' @param regn.cM disctance in cM form QTL for defining non-recombinanct mice.
##' @param max.p max number of expression traits used.
##'
add.geno <- function(transbands, cross, regn.cM=5, max.p=100){
  allele <- attr(cross, "alleles")
  gt <- c(paste0(allele[1], allele[1]),
          paste0(allele[1], allele[2]),
          paste0(allele[2], allele[2]))
  for(i in 1:length(transbands)){
    chr <- attr(transbands[[i]], "info")$chr
    Y <- attr(transbands[[i]], "Y")

    ## scan the first 100.
    Y <- Y[, 1:min(max.p, ncol(Y))]

    out.scan1 <- scanone.mvn(cross, Y, chr=chr)
    pos <- max(out.scan1)$pos
    pmarker <- find.pseudomarker(cross, chr=chr, pos=pos, where="prob")
    genoprob <- pull.genoprob(cross, chr=chr)
    genoprob.pmarker <- genoprob[, paste(pmarker, gt, sep=":")]
    geno <- apply(genoprob.pmarker, 1, which.max)

    m <- which(out.scan1$pos >= pos-regn.cM & out.scan1$pos <= pos+regn.cM)
    g <- apply(cross$geno[[chr]]$prob[,m,], 1:2, which.max)
    nonrecomb <- which(sapply(apply(g, 1, unique), length) == 1)
    names(nonrecomb) <- rownames(Y)[nonrecomb]

    attr(transbands[[i]], "geno") <- geno
    attr(transbands[[i]], "nonrecomb") <- nonrecomb
    attr(transbands[[i]], "info")$pos <- pos
    attr(transbands, "trans.info")$pos[i] <- pos
  }
  return(transbands)
}

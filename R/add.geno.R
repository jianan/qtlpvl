##' find the peak of scanone.mvn and the QTL genotype.
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

    pmarker1 <- find.pseudomarker(cross, chr=chr, pos=pos-regn.cM, where="prob")
    pmarker2 <- find.pseudomarker(cross, chr=chr, pos=pos+regn.cM, where="prob")
    geno1 <- apply(genoprob[, paste(pmarker1, gt, sep=":")], 1, which.max)
    geno2 <- apply(genoprob[, paste(pmarker2, gt, sep=":")], 1, which.max)
    nonrecomb <- which(geno==geno1 & geno==geno2)
    
    attr(transbands[[i]], "geno") <- geno
    attr(transbands[[i]], "nonrecomb") <- nonrecomb
    attr(transbands[[i]], "info")$pos <- pos
    attr(transbands, "trans.info")$pos[i] <- pos
  }
  return(transbands)
}

##' find the QTL addtive and dominant effect for all traits in the transband.
add.effect <- function(transbands, cross){
  allele <- attr(cross, "alleles")
  gt <- c(paste0(allele[1], allele[1]),
          paste0(allele[1], allele[2]),
          paste0(allele[2], allele[2]))
  for(i in 1:length(transbands)){  
    chr <- attr(transbands[[i]], "info")$chr
    Y <- attr(transbands[[i]], "Y")
    genoprob <- pull.genoprob(cross, chr=chr)
    out <- attr(transbands[[i]], "out")
    out$pmarker <- find.pseudomarker(cross, chr=chr, pos=out$pos, where="prob")

    p <- ncol(Y)
    eff1 <- eff2 <- eff3 <- numeric(p)
    for(j in 1:p){
      pmarker <- out[j, "pmarker"]
      genoprob.pmarker <- genoprob[, paste(pmarker, gt, sep=":")]
      geno <- apply(genoprob.pmarker, 1, which.max)
      eff1[j] <- mean(Y[geno==1, j], na.rm=TRUE)
      eff2[j] <- mean(Y[geno==2, j], na.rm=TRUE)
      eff3[j] <- mean(Y[geno==3, j], na.rm=TRUE)
    }
    out$eff.a <- (eff3 - eff1)/2
    out$eff.d <- eff2 - (eff3 + eff1)/2

    attr(transbands[[i]], "out") <- out[, c("pheno", "chr", "pos", "lod1", "eff.a", "eff.d")]
  }
  return(transbands)
}

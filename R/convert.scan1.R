##' Summarize result of scanone
##'
##' Turn the result of scanone to a smaller and thinner dataframe with
##' max LOD score and its position for each trait and each chrosome.
##'
##' @param out output from R/scanone
##' @param phenoname names of the traits of interest
##' @param chr chromosomes  of interest
##' @return data frame with \code{pheno, chr, pos}, and \code{lod1}, a
##' summary of scan1 result, max LOD score and its position for each
##' phenoname and chromosome combinantion.
##' @export
convert.scan1 <- function(out, phenoname, chr=1:19){
  CHR <- out[, "chr"]
  POS <- out[, "pos"]

  out <- out[, phenoname]
  maxLOD <- matrix(NA, p, length(chr))
  colnames(maxLOD) <- chr
  rownames(maxLOD) <- phenoname
  maxPOS <- maxLOD

  for(i in 1:length(chr)){
    index <- CHR==chr[i]
    LOD <- out[index, ]
    maxLOD[, i] <- apply(LOD, 2, max)
    maxPOS[, i] <- POS[index][apply(LOD, 2, which.max)]
  }

  out1 <- data.frame(pheno = rep(phenoname, ncol(maxLOD)),
                     chr = rep(colnames(maxLOD), each=nrow(maxLOD)),
                     lod1 = c(maxLOD),
                     stringsAsFactors = FALSE)
  out1$pos <- c(maxPOS)
  return(out1)
}

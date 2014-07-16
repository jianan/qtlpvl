##' For a scan1 result and probepos, find if the eQTLs are cis-eQTL or trans-eQTL.
##'
##' @param out1 summary of scan1 result
##' @param probepos Genetic position for each probe in cM.
##' @param chr Optional vector indicating the chromosomes considered.
##' @param marker.info Data.frame of marker information, result of
##' \code{get.marker.info}
##' @param lod.thr LOD threshold
##' @param trans.cM threshold of the distance between the qtl position
##' and the probe position to call a eQTL cis.
##' @return input data.frame with two added columns, "is.cis" and
##' "is.trans".
##' @export
get.trans.info <- function(out1, probepos, chr=NULL, marker.info, lod.thr=5, trans.cM=10){
  
  stopifnot(c("chr", "cM") %in% colnames(probepos))
  stopifnot(c("pheno", "chr", "pos", "lod1") %in% colnames(out1))

  if(!is.null(chr))  out <- subset(out1, chr %in% chr)
  out$chr0 <- probepos[out$pheno, "chr"]
  out$pos0 <- probepos[out$pheno, "cM"]
  out$is.cis <- (out$lod1 > lod.thr) & (out$chr0 == out$chr) & abs(out$pos0 - out$pos) < trans.cM
  out$is.trans <- (out$lod1 > lod.thr) & !out$is.cis
  
  return(out)
}

##' From a scan1 result, first count how many trans-eQTL are there in
##' every bin across the chromosomes, then find phenos inside each
##' transband. Add all pheno and geno infomations needed.
##'
##' @inheritParams find.trans
##' @inheritParams add.geno
##' @inheritParams scanone.mvn
##' @param mlratio matrix of expression traits
##' @return a list of transband, with each transband contains the
##' infomations needed for exploration plots.
##'
##' @export
make.transbands <- function(out1, probepos, cross, chr, marker.info, mlratio,
                            lod.thr=5, trans.cM=5, kernal.width=1, window.cM=10,
                            trans.count.thr=300, regn.cM=5){
  if(missing(marker.info))
      marker.info <- get.marker.info(cross, chr)
  transbands <- find.trans(out1, probepos, chr, marker.info,
                           lod.thr=lod.thr, trans.cM=trans.cM,
                           kernal.width=kernal.width, window.cM=window.cM,
                           trans.count.thr=trans.count.thr)
  transbands <- add.pheno(transbands, mlratio)
  transbands <- add.geno(transbands, cross, regn.cM=regn.cM)
  transbands <- add.effect(transbands, cross)
  return(transbands)
}

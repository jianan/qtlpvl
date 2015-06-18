##' From a scan1 result, count how many trans-eQTL are there in every
##' bin across the chromosomes.
##'
##' @inheritParams get.trans.info
##' @param kernal.width for weighted count, kernal.width is used as sd
##' value to generate normal weights. These weights are then used as
##' kernal function to sum over each bin/window.
##' @param window.cM bin width in cM.
##' @return a data.frame with "chr", "pos", "count", "count.sum" and
##' "count.wtsum".
##' @export
count.trans <- function(out1, probepos, chr, marker.info,
                        lod.thr=5, trans.cM=5, kernal.width=1, window.cM=10){

  stopifnot(c("pheno", "chr", "pos", "lod1") %in% colnames(out1))
  stopifnot(c("chr", "cM") %in% colnames(probepos))
  marker.info <-marker.info[marker.info$chr %in% chr, ]

  if(!"is.trans" %in% colnames(out1)){
    out1 <- get.trans.info(out1, probepos, chr, marker.info,
                           lod.thr=lod.thr, trans.cM=trans.cM)
  }
  out1 <- subset(out1, is.trans)
  if(nrow(out1) == 0) return()

  for(i.chr in chr){
    n.marker <- sum(marker.info$chr == i.chr)
    count <- numeric(n.marker)
    pos <- marker.info[marker.info$chr==i.chr, "pos"]
    for(j in 1:n.marker){
      count[j] <- sum(out1$is.trans & out1$chr == i.chr & out1$pos == pos[j])
    }
    count.sum <- numeric(n.marker)
    count.wtsum <- numeric(n.marker)
    for(j in 1:n.marker){
      dis <- abs(pos[j] - pos)
      index <- which(dis < window.cM/2)
      dis <- dis[index]
      num <- count[index]
      count.sum[j] <- sum(num)
      wt <- dnorm(dis, mean=0, sd=kernal.width)
      wt <- wt/dnorm(0, mean=0, sd=kernal.width)
      count.wtsum[j] <- sum(num*wt)
    }
    marker.info[marker.info$chr==i.chr, "count.sum"] <- count.sum
    marker.info[marker.info$chr==i.chr, "count.wtsum"] <- count.wtsum
    marker.info[marker.info$chr==i.chr, "count"] <- count
  }
  attr(marker.info, "class") <- c("scanone", "data.frame")

  return(marker.info)
}

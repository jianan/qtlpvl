##' From a scan1 result, first count how many trans-eQTL are there in
##' every bin across the chromosomes, then find phenos inside each
##' transband.
##' 
##' @inheritParams count.trans
##' @param trans.count.thr threshold used to define a region as
##' transband.
##' @return a list of transband.
##' @export
find.trans <- function(out1, probepos, chr, marker.info, 
                       lod.thr=5, trans.cM=5, kernal.width=1, window.cM=10,
                       trans.count.thr=300){
  
  stopifnot(c("pheno", "chr", "pos", "lod1") %in% colnames(out1))
  stopifnot(c("chr", "cM") %in% colnames(probepos))
  marker.info <- marker.info[marker.info$chr %in% chr, ]
  
  if(!"is.trans" %in% colnames(out1)){
    out1 <- get.trans.info(out1, probepos, chr, marker.info,
                           lod.thr=lod.thr, trans.cM=trans.cM)
  }
  out1 <- subset(out1, is.trans)
  if(nrow(out1) == 0) return()
  
  ## find trans.info, chr, left and right bound of each transband.
  out <- count.trans(out1, probepos, chr, marker.info, 
                     lod.thr=lod.thr, trans.cM=trans.cM,
                     kernal.width=kernal.width, window.cM=window.cM)
  out <- out[out$count.sum > trans.count.thr, ]

  ## x <- data.frame(pos=c(1:5, 15:20, 50:57),
  ##                 count=0)
  ## x$count <- sample(1:100, size=nrow(x))

  trans.info <- ddply(out, .(chr), function(x) {
    ## if two peaks are near each other(<5cM), combine them;
    ## otherwise treat them as different transband. 
    if(all(diff(x$pos) < trans.cM)) {
      res <- matrix(c(range(x$pos), x$pos[which.max(x$count)]), 1,3)
      return(res)
    }else{ 
      z <- which(diff(x$pos) >= trans.cM)
      z <- c(0, z, length(x$pos))
      nz <- length(z)
      res <- matrix(NA, nrow=(nz-1), ncol=3)
      for(i in 1:(nz-1)){
        pos <- x$pos[(1+z[i]):z[i+1]]
        count <- x$count[(1+z[i]):z[i+1]]
        res[i, ] <- c(range(pos), pos[which.max(count)])
      }
      return(res)
    }})
  ## peak is the position with the most trans-eQTLs
  colnames(trans.info) <- c("chr", "trans.l", "trans.r", "peak")

  ## expand the transband, +- 5cM on the ranges
  trans.info$trans.l <- trans.info$trans.l - window.cM/2
  trans.info$trans.l[trans.info$trans.l < 0] <- 0
  trans.info$trans.r <- trans.info$trans.r + window.cM/2
  if(nrow(trans.info) == 0) return()

  ## find phenos inside each transband.
  trans.pheno <- list()
  for(i in 1:nrow(trans.info)){
    chr <- trans.info[i, "chr"]
    trans.l <- trans.info[i, "trans.l"]
    trans.r <- trans.info[i, "trans.r"]
    peak <- trans.info[i, "peak"]
    out.trans <- out1[out1$chr==chr & out1$pos > trans.l & out1$pos < trans.r, ]
    out.trans <- out.trans[order(out.trans$lod1, decreasing=TRUE), ]
    trans.pheno[[i]] <- out.trans$pheno
    trans.info$count[i] <- nrow(out.trans)
    attr(trans.pheno[[i]], "info") <- list(chr=chr, trans.l=trans.l, trans.r=trans.r,
                                           peak=peak, count=nrow(out.trans))
    attr(trans.pheno[[i]], "out") <- out.trans
  }
  rownames(trans.info) <- names(trans.pheno) <-
      paste0("chr", trans.info$chr, ".", as.integer(trans.info$peak))
  attr(trans.pheno, "trans.info") <- trans.info
  return(trans.pheno)
}

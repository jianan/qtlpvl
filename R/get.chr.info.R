get.chr.info <- function(marker.info){
  chr <- unique(marker.info$chr)
  start <- NULL
  end <- NULL
  for(i in chr){
    start[i] <- min(marker.info[marker.info$chr==i, "pos"])
    end[i] <- max(marker.info[marker.info$chr==i, "pos"])
  }
  len <- end-start
  data.frame(chr, start, end, len, stringsAsFactors=FALSE)
}

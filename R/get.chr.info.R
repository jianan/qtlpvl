##' @export
get.chr.info <- function(marker.info){
  require(plyr) || stop(" the required package 'plyr' is not installed. ")
  chr <- unique(marker.info$chr)
  chr.info <- ddply(marker.info, .(chr), summarise,
                    start=min(pos),
                    end=max(pos))
  rownames(chr.info) <- chr.info$chr
  chr.info <- chr.info[chr, ]
  start <- chr.info$start
  end <- chr.info$end
  len <- end-start
  result <- data.frame(chr, start, end, len, stringsAsFactors=FALSE)
  return(result)
}

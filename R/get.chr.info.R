##' @importFrom plyr ddply summarise
get.chr.info <- function(marker.info){
  chr <- unique(marker.info$chr)
  chr.info <- ddply(marker.info, .(chr), summarise,
                    start=min(pos),
                    end=max(pos))
  rownames(chr.info) <- chr.info$chr
  chr.info <- chr.info[chr, ]
  start <- chr.info$start
  end <- chr.info$end
  len <- end-start
  data.frame(chr, start, end, len, stringsAsFactors=FALSE)
}

# avoid warnings from R CMD check
globalVariables(c("pos", "."))

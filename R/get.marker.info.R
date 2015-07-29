##' Marker information.
##'
##' Pull out marker names positions and chromosomes as a data frame.
##'
##' @param cross An object of class \code{cross}. See
##' \code{read.cross} for details.
##' @param chr Optional vector indicating the chromosomes for which
##' markers info need to be pulled out. This should be a vector of
##' character strings referring to chromosomes by name; numeric values
##' are converted to strings.
##' @return A data.frame whose first column contains the chromosome
##' IDs, second column contains cM positions, row names are
##' markernames.
##'
##' @examples
##' data(fake.phenos)
##' marker <- get.marker.info(listeria)
##'
##' @export
get.marker.info <- function(cross, chr){
  if (!any(class(cross) == "cross"))
      stop("Input should have class \"cross\".")
  if (!missing(chr))
      cross <- subset(cross, chr = chr)
  out <- NULL
  for(i.chr in names(cross$geno)){
    map.chr <- attr(cross[[c("geno", i.chr, "prob")]], "map")
    ind <- substr(names(map.chr),1,3) =="loc"
    names(map.chr)[ind] <- paste0("c", i.chr, ".",names(map.chr)[ind])
    out <- rbind(out, data.frame(chr=i.chr, pos=map.chr, stringsAsFactors=FALSE))
  }
  return(out)
}

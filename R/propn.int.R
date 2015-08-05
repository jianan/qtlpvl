##' propotion plot.
##'
##' Compare pred.test to the genotype of each marker, calculate and
##' plot thr propotion of matches at each marker position and inside
##' each marker interval, could be used to estimate interval for the
##' common eQTL in a transband.
##'
##' @inheritParams mosaic.plot
##' @param thr cutoff for the interval. default is 0.9
##' @return infered interval.
##' @export
propn.int <- function(genotype, map, pred.test, pred.score, weighted = FALSE, thr=0.9){

  stopifnot(is.matrix(genotype))
  genotype.match <- genotype == pred.test
  n <- length(pred.test)
  m <- length(map)
  stopifnot(nrow(genotype) == n & ncol(genotype) == m)

  if(!weighted){
    propn <- apply(genotype.match, 2, mean, na.rm=TRUE)
  }else{
    propn <- apply(genotype.match * pred.score, 2, mean, na.rm=TRUE)
  }

  ind <- which(propn > thr)
  if(length(ind)==0) return()
  int <- c("int.l" = map[max(1, min(ind)-1)], "int.r" = map[min(m, max(ind)+1)])
  return(int)
}

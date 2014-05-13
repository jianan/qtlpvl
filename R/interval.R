##' find interval for multiple traits mapped to the same region.
##'
##' @inheritParams group.train.test
##' @inheritParams mosaic.plot
##' @param do.plot logical, plot results? default is \code{TRUE}
##' @return list of chr and interval and plots if \code{do.plot=TRUE}.
##' @export
interval <- function(cross, Y, chr, region.l, region.r, do.plot=TRUE, main="", weighted=FALSE, thr=0.9){
  
  group <- group.train.test(cross, Y, chr, region.l, region.r)
  data.train <- group$data.train
  data.test <- group$data.test
  class.train <- group$geno.train
  genotype <- group$geno.test
  map <- group$map

  fit <- classification(data.train, data.test, class.train, method="LDA")
  pred.test <- fit$pred.test
  pred.score <- fit$pred.score
  sca <- fit$sca
  error.train <- fit$error.train

  int <- propn.int(genotype, map, pred.test, pred.score, weighted = weighted, thr=thr)
  if(do.plot){
    layout(matrix(c(1,2,3,3),2,2))
    mosaic.plot(genotype, map, pred.test, pred.score, main=main, weighted = weighted, label=FALSE)
    propn.plot(genotype, map, pred.test, pred.score, weighted = weighted)
    plotlda(data.train, data.test, class.train, pred.test, error.train, sca)
  }
  
  result <- list(chr=chr, int=int, map=map)
  return(result)
}

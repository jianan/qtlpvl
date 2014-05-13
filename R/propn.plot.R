##' propotion plot.
##'
##' Compare pred.test to the genotype of each marker, calculate and
##' plot the propotion of matches at each marker position and inside
##' each marker interval, could be used to estimate interval for the
##' common eQTL in a transband.
##' 
##' @inheritParams mosaic.plot
##' @return a proportion plot.
##' @export
propn.plot <- function(genotype, map, pred.test, pred.score, main="", weighted = FALSE){
  
  stopifnot(is.matrix(genotype))
  genotype.match <- genotype == pred.test
  n <- length(pred.test)
  m <- length(map)
  stopifnot(nrow(genotype) == n & ncol(genotype) == m)
  
  if(!weighted){
    propn <- apply(genotype.match, 2, mean, na.rm=TRUE)
    ylab <- "proportion"
  }else{
    propn <- apply(genotype.match * pred.score, 2, mean, na.rm=TRUE)
    ylab <- "weighted proportion"
  }
  
  plot(map, propn, type="b", ylim=c(0, 1), pch=20, 
       xaxt="n", xlab="marker position (cM)", 
       ylab=ylab, main=main)
  axis(1, at=map, labels = sprintf("%1.1f", map))
  abline(h=1, lty=3)
  abline(v=map, lty=3)
  
  ## insert a point between markers, only count when left==right.
  map.ins <- numeric(m-1)
  genotype.ins <- matrix(NA, nrow=n, ncol=m-1)
  for(i in 1:(m-1)){
    map.ins[i] <- (map[i] + map[i+1])/2
    for(j in 1:n){
      if(genotype.match[j, i] == genotype.match[j, i+1]){
        genotype.ins[j, i] <- genotype.match[j, i]
      }
    }
  }

  if(!weighted){
    propn <- apply(genotype.ins, 2, mean, na.rm=TRUE)
  }else{
    propn <- apply(genotype.ins * pred.score, 2, mean, na.rm=TRUE)
  }

  points(map.ins, propn, type="p", pch=20, col="red")
} 

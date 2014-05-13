##' mosaic plot
##'
##' Comparasion of pred.test to the genotype of each marker, used to
##' estimate interval for the common eQTL in a transband.
##'
##' @param genotype genotype matrix, rows are samples, columns are markers.
##' @param map genetic map for markers.
##' @param pred.test vector of length n, predicted class (1, 2 or 3)
##' for each sample.
##' @param pred.score measure of how good is the prediction, when
##' weighted = TRUE, this score is used to order the samples.
##' @param main main titile of the plot.
##' @param weighted use pred.score as weights? default is \code{FALSE}.
##' @param label include sample name as labels in the plot? default is
##' \code{TRUE}.
##' @return a mosaic plot.
##' @export
mosaic.plot <- function(genotype, map, pred.test, pred.score, main="", weighted = FALSE, label=TRUE){

  stopifnot(is.matrix(genotype))
  genotype.match <- genotype == pred.test
  n <- length(pred.test)
  m <- length(map)
  stopifnot(nrow(genotype) == n & ncol(genotype) == m)
  
  ## for each sample, region of exact match.
  region.l <- apply(genotype.match, 1, function(x){suppressWarnings(min(which(x)))})
  region.r <- apply(genotype.match, 1, function(x){suppressWarnings(max(which(x)))})

  if(!weighted){
    o <- order(region.r, 0-region.l, decreasing=FALSE)
  }else{
    brks <- seq(from = 0.5, to=1, by=0.1)
    pred.brks <- as.numeric(cut(pred.score, breaks=c(0,brks)))
    n.brks <- 0
    for(i in 1:length(brks)){
      n.brks[i] <- sum(pred.brks <= i)
    }
    o <- order(pred.brks, region.r, 0-region.l, decreasing=FALSE)
  }
  
  genotype.match <- genotype.match[o, ]
  
  plot(0, 0, xlim=range(map), ylim=c(1, n + 1), 
       xaxt="n", type="n", xlab="marker", ylab="Mouse", main=main)
  axis(1, at=map, labels = sprintf("%1.1f", map))

  if(label){
    mouse <- rownames(genotype.match)
    axis(side=4, at=1:n, labels = mouse, tick=FALSE, las=1)
  }
  
  for(i in 1:n){
    rect(xleft=map[1], xright=map[m], ybottom=i, ytop=i+0.6, col="grey")
    for(marker.i in 1:(m-1)){
      if(genotype.match[i, marker.i] == genotype.match[i, marker.i+1]){ 
        rect(xleft=map[marker.i], xright=map[marker.i +1], ybottom=i, ytop=i+0.6, 
             col=c("red", "green")[1+genotype.match[i, marker.i]])
      }
    }
  }
  
  if(weighted){
    abline(h=n.brks[1:(length(brks)-1)]+1, col="red", lty=2)
    for(i in 1:(length(brks)-1)){
      axis(side=4, at=1+n.brks[i], labels=brks[i+1], col="red", cex=0.7)
    }
  }
  
}

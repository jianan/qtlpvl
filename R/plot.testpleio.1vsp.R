##' @export
plot.testpleio.1vsp <- function(x, ...){
  if (!any(class(x) == "testpleio.1vsp")) 
      stop("Input should have class \"testpleio.1vsp\".")
  dots <- list(...)
  
  LODdiff <- x$LODdiff
  LOD1 <- x$LOD1
  map <- x$map
  map.marker <- x$map.marker
  maxLOD <- x$maxLOD
  maxPOS <- x$maxPOS
  rg <- range(map)
  
  plot(y=LOD1, x=map, type="l",
       ylim=c(0,max(LOD1,na.rm=TRUE)), xlim=rg,
       xaxt="n", xlab="cM position", ylab="LOD")
  rug(map.marker, ticksize=-0.01)
  axis(side=1, at=map.marker, labels=sprintf("%.1f",map.marker))
  points(maxPOS, maxLOD, pch=20) ## single trait result.
  points(x=map[which.max(LOD1)], y=max(LOD1),
         col="blue", pch=17) ## plot the best one QTL
}

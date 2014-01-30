##' @export
plot.testpleio.1vs2 <- function(x, ...){
  
  object <- x
  if (!any(class(object) == "testpleio.1vs2")) 
      stop("Input should have class \"testpleio.1vs2\".")
  dots <- list(...)

  LOD1 <- object$LOD1 ## LOD of joint mapping 
  LOD2 <- object$LOD2 ## 2 QTL, LOD score matrix.
  LODdiff <- object$LODdiff
  maxPOS <- object$maxPOS ## single trait mapped pos
  maxLOD <- object$maxLOD ## single trait LOD
  Group <- object$Group
  pvalue <- object$pvalue
  map <- object$map
  map.marker <- object$map.marker
  rg <- range(map)

  old.yaxt <- par("yaxt")
  old.mfrow <- par("mfrow")
  on.exit(par(yaxt = old.yaxt, mfrow = old.mfrow))
  
  par(mfrow=c(2,1))
  plot(y=LOD1, x=map, type="l",
       ylim=c(0,max(LOD2,na.rm=TRUE)), xlim=rg,
       xaxt="n", xlab="cM position", ylab="LOD")
  rug(map.marker, ticksize=-0.01)
  axis(side=1, at=map.marker, labels=sprintf("%.1f",map.marker))
  points(maxPOS, maxLOD, col=c("blue","red")[Group], pch=20) ## single trait result.
  points(x=attr(LODdiff, "LOD1pos"), y=attr(LODdiff, "LOD1lod"),
         col="black", pch=17) ## plot the best one QTL
  points(y=c(0,0), x=attr(LODdiff, "LOD2pos"),
         col=c("blue","red"), pch=2)  ## plot the best two QTLs

  ind <- arrayInd(which.max(LOD2), .dim=dim(LOD2))
  points(map,LOD2[, ind[2]],type="l",col="blue")
  points(map,LOD2[ind[1], ],type="l",col="red")
  
  plot(x=1:(length(Group)-1), y=object$LODdiff.trace,
       xlab="i.cut", ylab="LODdiff", type="b")
  
}

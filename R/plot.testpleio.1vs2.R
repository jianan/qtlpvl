##' @export
plot.testpleio.1vs2 <- function(x, xlab="Map position (cM)", ylab="LOD",
                                mgp=c(1.6, 0.2, 0),
                                ...){
  
  if (!any(class(x) == "testpleio.1vs2")) 
      stop("Input should have class \"testpleio.1vs2\".")

  LOD1 <- x$LOD1 ## LOD of joint mapping 
  LOD2 <- x$LOD2 ## 2 QTL, LOD score matrix.
  LODdiff <- x$LODdiff
  maxPOS <- x$maxPOS ## single trait mapped pos
  maxLOD <- x$maxLOD ## single trait LOD
  Group <- x$Group
  pvalue <- x$pvalue
  map <- x$map
  map.marker <- x$map.marker
  rg <- range(map)
  ylim <- c(min(0, min(maxLOD, LOD1, LOD2, na.rm=TRUE)),
            max(maxLOD, LOD1, LOD2, na.rm=TRUE))

  plot(y=LOD1, x=map, type="l",
       ylim=ylim, xlim=rg, xaxt="n", xlab=xlab, ylab=ylab,
       mgp=mgp, ...)

  rug(map.marker, ticksize=-0.01)
  axis(side=1, at=map.marker, labels=sprintf("%.1f",map.marker), mgp=mgp)
  points(maxPOS, maxLOD, col=c("blue","red")[Group], pch=20) ## single trait result.
  points(x=attr(LODdiff, "LOD1pos"), y=attr(LODdiff, "LOD1lod"),
         col="black", pch=17) ## plot the best one QTL
  points(y=c(0,0), x=attr(LODdiff, "LOD2pos"),
         col=c("blue","red"), pch=2)  ## plot the best two QTLs

  ## add LOD profile for the two QTLs.
  ind <- arrayInd(which.max(LOD2), .dim=dim(LOD2))
  points(map,LOD2[, ind[2]],type="l",col="blue")
  points(map,LOD2[ind[1], ],type="l",col="red")
}

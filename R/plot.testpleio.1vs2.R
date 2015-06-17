##' @export
plot.testpleio.1vs2 <- function(x, xlab="Map position (cM)", ylab="LOD",
                                mgp=c(1.6, 0.2, 0), ...){

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
  ylim <- c(0, max(maxLOD, LOD1, LOD2, na.rm = TRUE)*1.05)
  broman::grayplot(y = LOD1, x = map, type = "l", ylim = ylim, xlim = rg,
                   xaxt = "n", xlab = xlab, ylab = ylab, mgp = mgp, lwd=2,
                   yat=pretty(ylim), yaxs="i", ...)
  rug(map.marker, ticksize = -0.01)
  points(maxPOS, maxLOD, col = c("slateblue", "violetred")[Group], pch = 20) ## single trait result.
  points(x = attr(LODdiff, "LOD1pos"), y = attr(LODdiff, "LOD1lod"),
         col = "black", pch = 17) ## plot the best one QTL
  points(y = c(2, 2), x = attr(LODdiff, "LOD2pos"),
         col = c("slateblue", "violetred"), pch = 17) ## add LOD profile for the two QTLs.
  ind <- arrayInd(which.max(LOD2), .dim = dim(LOD2))
  lines(map, LOD2[, ind[2]], col = "slateblue", lwd=2)
  lines(map, LOD2[ind[1], ], col = "violetred", lwd=2)
}

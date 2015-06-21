##' @export
plot.testpleio.1vsp <- function(x, xlab="Map position (cM)", ylab="LOD score",
                                mgp=c(1.6, 0.2, 0), ...){

  if (!any(class(x) == "testpleio.1vsp"))
      stop("Input should have class \"testpleio.1vsp\".")

  LODdiff <- x$LODdiff
  LOD1 <- x$LOD1
  map <- x$map
  map.marker <- x$map.marker
  maxLOD <- x$maxLOD
  maxPOS <- x$maxPOS
  rg <- range(map)

  ylim <- c(0, max(maxLOD, LOD1, na.rm = TRUE)*1.05)
  broman::grayplot(y = LOD1, x = map, type = "l", ylim = ylim, xlim = rg,
                   xaxt = "n", xlab = xlab, ylab = ylab, mgp = mgp, lwd=2,
                   yat=pretty(ylim), yaxs="i", ...)
  rug(map.marker, ticksize = -0.01)
  points(maxPOS, maxLOD, pch = 20) ## single trait result.
  points(x=map[which.max(LOD1)], y=max(LOD1),
         col="slateblue", pch=17) ## plot the best one QTL
}

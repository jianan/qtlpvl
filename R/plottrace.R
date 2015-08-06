##' plot trace
##'
##' plot trace of a searching method to see how the object function
##' changes as the search goes.
##'
##' @param x object to plot
##' @param ... parameters to be passed through to plotting functions.
##' @export
plottrace <- function(x, ...)
    UseMethod("plottrace")

##' plot trace
##'
##' plot trace of \code{testpleio.1vs2} to see how \code{LODdiff}
##' changes as the cutting point moves.
##'
##' @param x object of class "testpleio.1vs2"
##' @param ... parameters to be passed through to plotting functions.
##' @param xlab,ylab,mgp same as in \code{plot}
##' @export
plottrace.testpleio.1vs2 <- function(x, xlab="cut point", ylab=expression(LOD["2v1"]),
                        mgp=c(1.6, 0.2, 0), ...){
  object <- x
  if (!any(class(object) == "testpleio.1vs2"))
      stop("Input should have class \"testpleio.1vs2\".")
  t <- object$LODdiff.trace
  t[t<0] <- 0
  n.max <- which.max(t)
  v.max <- max(t, na.rm=TRUE)
  ylim <- c(-v.max*0.05, v.max*1.05)
  broman::grayplot(1:length(t), t, type = "n", xlab = xlab, ylab = ylab,
                   ylim = ylim, yaxs="i", mgp = mgp,
                   xlim=c(0.5, length(t)+0.5), xaxs="i", las=1, ...)
  u <- par()$usr
  segments(x0=n.max, y0=u[3], x1=n.max, y1=v.max, lwd=2, col="violetred")
  lines(seq(along=t), t, type="o", pch=21, bg="slateblue", cex=0.8)
}

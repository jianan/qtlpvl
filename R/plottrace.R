##' plot trace
##'
##' plot trace of a searching method to see how the object function
##' changes as the search goes.
##'
##' @param x object to plot
##' @export
plottrace <- function(x)
    UseMethod("plottrace")

##' plot trace
##'
##' plot trace of \code{testpleio.1vs2} to see how \code{LODdiff}
##' changes as the cutting point moves.
##'
##' @param x object of class "testpleio.1vs2"
##' @export
plottrace.testpleio.1vs2 <- function(x){
  
  object <- x
  if (!any(class(object) == "testpleio.1vs2")) 
      stop("Input should have class \"testpleio.1vs2\".")

  plot(x=1:(length(object$Group)-1), y=object$LODdiff.trace,
       xlab="i.cut", ylab="LODdiff", type="b")
  n.max <- which.max(object$LODdiff.trace)
  v.max <- max(object$LODdiff.trace, na.rm=TRUE)
  u <- par()$usr
  segments(x0=n.max, y0=u[3], x1=n.max, y1=v.max, lty=3, col="red")
}

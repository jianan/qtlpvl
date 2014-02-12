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
}

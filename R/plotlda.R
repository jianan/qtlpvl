##' plot of ld2 vs ld1.
##'
##' plot training and testing data together, to compare their
##' pattern. different pattern indicates multiple QTLs, because it
##' needs extra signal to drag the tesing points away from the traing
##' points.
##'
##' @inheritParams mosaic.plot
##' @inheritParams classification
##' @param error.train LDA traning error
##' @param sca scalling factor for ld1 and ld2
##' @return plot of ld2 vs ld1
##' @export
plotlda <- function(data.train, data.test, class.train,
                    pred.test, error.train, sca, main=""){

  n.train <- nrow(data.train)
  n.test <- nrow(data.test)
  disc.test <- as.matrix(data.test) %*% sca
  disc.train <- as.matrix(data.train) %*% sca
  xlim <- range(disc.test[,1], disc.train[,1])
  ylim <- range(disc.test[,2], disc.train[,2])
  ## ylim <- xlim <- range(disc.test, disc.train)
  
  plot(disc.test, col=c(2, 3, 4)[pred.test], pch=c("1", "2", "3")[pred.test], 
       cex=0.7, xlim=xlim, ylim=ylim, main=main, xlab="LD1", ylab="LD2")
  text <- paste0("non-recomb=", n.train, ", recomb=", n.test, ## sample size
                 ", error.train=", round(error.train, digits=3))
  mtext(text, side=3, line=0)
  points(disc.train, col="grey", pch=c("1", "2", "3")[class.train], cex=0.7)
  
}

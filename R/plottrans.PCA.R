##' PCA plot for a transband.
##' 
##' @param Y matrix of traits in the transband
##' @param geno genotype at the peak marker of scanone.mvn. 
##' @param nonrecomb index for non-recombnent(in a 10cM region) individuals.
##' @param max.p max number of traits used for PCA analysis. This is
##' used to avoid rank deficiency.
##' @export
plottrans.PCA <- function(Y, geno, nonrecomb, max.p=100, ...){

  ## use the first 100.
  Y <- Y[, 1:min(max.p, ncol(Y))]

  PC <- predict(princomp(Y))[, 1:2]
  
  plot(PC, xlab="PC1", ylab="PC2", ...)
  points(PC[-nonrecomb, ], pch=20, col="yellow")
  points(PC[nonrecomb, ], pch=20, col=geno[nonrecomb])
}

##' LDA plot for a transband.
##'
##' @inheritParams plottrans.PCA
##' @export
plottrans.LDA <- function(Y, geno, nonrecomb, max.p=100, ...){
  suppressMessages(require(MASS))|| stop("the required package 'MASS' is not installed. ")

  ## use the first 100.
  Y <- Y[, 1:min(max.p, ncol(Y))]

  fit <- lda(Y[nonrecomb, ], grouping=geno[nonrecomb])
  pred <- predict(fit, Y[-nonrecomb, ])
  sca <- fit$scaling
  y.lda <- Y %*% sca

  plot(y.lda, col=geno, xlab="LD1", ylab="LD2", ...)
  points(y.lda[nonrecomb, ], pch=20, col=geno[nonrecomb])
  points(y.lda[-nonrecomb, ], pch=20, col="yellow")
}

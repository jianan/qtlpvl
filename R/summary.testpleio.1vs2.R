##' @export
print.summary.testpleio.1vs2 <- function(x, ...){
  
  LODdiff <- x$LODdiff
  Group <- x$Group
  cat("Single QTL model: ")
  cat("lod ")
  cat(sprintf("%.2f",attr(LODdiff,"LOD1lod")))
  cat(", pos ")
  cat(paste(sprintf("%.2f",attr(LODdiff,"LOD1pos"))," cM",sep=""))
  cat("\n")

  cat("Two QTLs model: ")
  cat("lod ")
  cat(sprintf("%.2f",attr(LODdiff,"LOD2lod")))
  cat(", pos ")
  cat(paste(paste(sprintf("%.2f",attr(LODdiff,"LOD2pos"))
                  ," cM",sep=""),collapse=", ", sep=""))
  cat("\n")
  cat("  Traits influenced by the left QTL: ")
  cat(which(Group==1))
  cat("\n")
  cat("  Traits influenced by the right QTL: ")
  cat(which(Group==2))
  cat("\n")
  
  cat("Difference of LOD score: ")
  cat(sprintf("%.2f",x$LODdiff))
  cat("\n")

  cat("P-value is: ")
  cat(x$pvalue)
  cat(" (from", x$n.simu ,"simulations)")
  cat("\n")
}


##' @export
summary.testpleio.1vs2 <- function(object, ...){
  if(!any(class(object) == "testpleio.1vs2"))
      stop("Input should have class 'testpleio.1vs2'. ")
  result <- object[c("LODdiff","Group","pvalue","n.simu")]
  class(result) <- "summary.testpleio.1vs2"
  return(result)
}

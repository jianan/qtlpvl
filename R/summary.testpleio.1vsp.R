##' @export
print.summary.testpleio.1vsp <- function(x, ...){
  
  LODdiff <- x$LODdiff
  LOD1 <- x$LOD1
  map <- x$map
  maxLOD <- x$maxLOD
  maxPOS <- x$maxPOS
  cat("Single QTL model: ")
  cat("lod ")
  cat(sprintf("%.2f",max(LOD1)))
  cat(", pos ")
  cat(paste(sprintf("%.2f",map[which.max(LOD1)])," cM",sep=""))
  cat("\n")

  cat("Multiple QTLs model: ")
  cat("lod ")
  cat(sprintf("%.2f",x$LODp))
  cat("\n")
  print(data.frame(Triat=1:length(maxLOD), POS=maxPOS, LOD=maxLOD))
  cat("\n")
  
  cat("Difference of LOD score: ")
  cat(sprintf("%.2f",x$LODdiff))
  cat("\n")

  cat("P-value is: ")
  cat(x$pvalue)
  cat(", (form", x$n.simu ,"simulations)")
  cat("\n")
}

##' @export
summary.testpleio.1vsp <- function(object, ...){
  if(!any(class(object) == "testpleio.1vsp"))
      stop("Input should have class 'testpleio.1vsp'. ")
  result <- object
  class(result) <- "summary.testpleio.1vsp"
  return(result)
}

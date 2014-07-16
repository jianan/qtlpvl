##' for each transband, add the phenotypes "Y".
add.pheno <- function(transbands, mlratio){
  for(i in 1:length(transbands)){
    phenos <- transbands[[i]]
    Y <- mlratio[, phenos]
    attr(transbands[[i]], "Y") <- Y
  }
  return(transbands)
}

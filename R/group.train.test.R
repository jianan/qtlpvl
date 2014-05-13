##' find the training and tesing set from genotypes in a given chr region.
##'
##' @param cross An object of class \code{cross}. See
##' \code{read.cross} for details.
##' @param Y traits mapped to the region.
##' @param chr which chromosome is under invest.
##' @param region.l left bound of region under invest.
##' @param region.r right bound of region under invest.
##' @return a list of data.train, data.test, geno.train, geno.test, map.
##' @export
group.train.test <- function(cross, Y, chr, region.l, region.r){

  if(!is.null(rownames(Y))){
    ID <- rownames(Y)
  }else if("id" %in% colnames(cross$pheno)){
    ID <- as.character(cross$pheno$id)
  }else{
    ID <- as.character(1:nrow(Y))
  }
  rownames(Y) <- ID
  
  if(length(chr) >1){
    chr <- chr[1]
    warning("Using only chr[1]")
  }
  cross <- cross[chr, ]
  
  ## fill in genotype
  cross.fill <- fill.geno(cross, method="argmax", error.prob=0.002, map.function="c-f")
  genotype <- pull.geno(cross.fill)
  rownames(genotype) <- ID
  
  map <- pull.map(cross)[[1]]
  genotype <- genotype[, map > region.l & map < region.r]
  if(ncol(genotype) < 2){
    stop("Need at least two markers inside the region.")
  }else{
    n.geno <- apply(genotype, 1, function(x){length(unique(x))})
    map <- map[map > region.l & map < region.r]
  }
  
  trainID <- ID[n.geno == 1]
  testID <- ID[n.geno != 1]
  
  data.train <- data.frame(Y[trainID, ])
  data.test <- data.frame(Y[testID, ])
  colnames(data.train) <- c(colnames(Y))
  colnames(data.test) <- colnames(Y)
  geno.train <- genotype[trainID, 1]
  geno.test <- genotype[testID, ]
  
  result <- list(data.train=data.train, data.test=data.test,
                 geno.train=geno.train, geno.test=geno.test, map=map)
  
  return(result)
}

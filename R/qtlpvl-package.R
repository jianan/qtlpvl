##' qtlpvl
##'
##' QTL: test of pleiotrophy vs. close linkage
##'
##' @import qtl
##' @import Rcpp
##' @import plyr
##' @name qtlpvl
##' @useDynLib qtlpvl
##' @docType package
NULL


##' Simulated Small Gene Expression Data Set
##'
##' Simulated gene expression data, with 10 expression traits for 120
##' individuals.
##'
##' @format A matrix with 120 individuals in rows, 10 traits in columns.
##'
##' @details These data were simulated using the genotype data from
##' the `listeria` data set (an F$_2$ population) provided with
##' R/qtl. The phenotypes were simulated using two markers from
##' chromosome 1 as QTL, with the first QTL having an additive allelic
##' effect, and with one of the alleles at the second QTL being
##' strictly dominant. There are 10 phenotypes. The first 5 are
##' controlled by the first QTL, and the other 5 traits are controlled
##' by the second QTL (and with a negative and larger effect). The 10
##' phenotypes were generated with these QTL effects plus independent,
##' normally distributed residual variation. Treating these traits as
##' gene expression measurements, we assigned genomic positions at
##' random. The phenotype data is stored in matrix `fake.phenos` and their
##' positions are stored in data frame `fake.probepos`.
##'
##' @docType data
##'
##' @keywords datasets
##'
##' @name fake.phenos
NULL


##' Simulated Genetic Positions for Gene Expression Data Set
##'
##' Simulated position for the 10 gene expression traits with
##' chromosome number and genetic positions in cM.
##'
##' @format A data frame with 10 rows and 2 columns: 'chr' and 'cM'.
##'
##' @details These data were simulated using the genotype data from
##' the `listeria` data set (an F$_2$ population) provided with
##' R/qtl. The phenotypes were simulated using two markers from
##' chromosome 1 as QTL, with the first QTL having an additive allelic
##' effect, and with one of the alleles at the second QTL being
##' strictly dominant. There are 10 phenotypes. The first 5 are
##' controlled by the first QTL, and the other 5 traits are controlled
##' by the second QTL (and with a negative and larger effect). The 10
##' phenotypes were generated with these QTL effects plus independent,
##' normally distributed residual variation. Treating these traits as
##' gene expression measurements, we assigned genomic positions at
##' random. The phenotype data is stored in matrix `fake.phenos` and their
##' positions are stored in data frame `fake.probepos`.
##'
##' @docType data
##'
##' @keywords datasets
##'
##' @name fake.probepos
NULL

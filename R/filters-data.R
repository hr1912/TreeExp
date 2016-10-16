
#' @title Filters for subsetting orthologous genes
#' 
#' @name genefilters
#' 
#' @description Two example filters for subsetting orthologous genes. 
#' 1. Housekeeping genes
#' 2. Nervous System Development genes 
#' 
#' @docType data
#' 
#' @format Two dataframes whose first column corresponds to Ensembl gene ids,
#' second column corresponds to housekeeping gene symbol,
#' and third column corresponds to Go ids related to Nervous System Development.
#' 
#' @references 
#' Brawand,D. et al. (2011) The evolution of gene expression levels in mammalian organs.
#' Nature, 478, 343-348.
#' Dorus S. et al. (2004) Accelerated evolution of nervous system genes in the origin of Homo sapiens. 
#' Cell 119:1027-1040.
#' Eisenberg E, Levanon EY. (2013) Human housekeeping genes, revisited. 
#' Trends Genet 29:569-574.
#' 
#' @examples 
#' data(genefilters)
#' head(tetrafilter)
#' 
NULL

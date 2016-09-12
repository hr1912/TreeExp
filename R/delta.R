#' @title Estimating delta values from a distance matrix.
#'
#' @name estdelta
#'
#' @rdname delta
#'
#' @param disMat a distance matrix with column and row names
#'
#' @return returns a list of two,
#' quartet: vector of delta numbers,
#' each corresponded to a quartet in the distance matrix.
#' taxa: matrix of delta numbers,
#' colnames corresponded to taxaNames
#' taxaNames: vector of taxaNames
#'
#' @details
#'
#' @examples
#' data(tetraexp)
#' dis.mat <- expdist(tetraexp.objects, taxa = "all",
#'                      subtaxa = "Brain",
#'                      method = "pea")
#' deltas <- estdelta(dis.mat)
#' hist(deltas$quartet)
#'
#' @references
#' Holland,B.R. et al. 2002. Delta plots: a tool for analyzing phylogenetic distance data.
#' Mol. Biol. Evol., 19, 2051-2059.
#'
#' @export
estdelta = function(disMat = NULL) {

  taxaNames <- row.names(disMat)
  taxaNumber <- length(taxaNames)

  allCombn <- combn(taxaNumber, 4)

  # for all the combinations  of quartet c(n,4)
  delta.quartet <- vector(mode = "numeric", length = length(allCombn[1,]))


  delta.taxa <- matrix(-1, nr = choose(taxaNumber-1,3), nc = taxaNumber)
  #names(delta.taxa) <- taxaNames

  #delta.quartet <- numeric(length=length(allCombn[1,]))
  for (i in 1:length(allCombn[1,])) {

    aCombn <- allCombn[,i]

    aQuartet <- disMat[aCombn,aCombn]
    aQuartet.melted <- melt(aQuartet)

    aQuartet.melted <- aQuartet.melted[aQuartet.melted$value!=0,]
    sixValues <- aQuartet.melted$value

    threeQuantities <- sort(c(sixValues[3]+sixValues[4],
                              sixValues[2]+sixValues[5],
                              sixValues[1]+sixValues[6]))

    delta.quartet[i] <- (threeQuantities[3] - threeQuantities[2]) / (threeQuantities[3] - threeQuantities[1])

    for (j in 1:4) {

      for (k in 1:length(delta.taxa[,aCombn[j]])) {

        if (delta.taxa[k, aCombn[j]] == -1) {
          delta.taxa[k, aCombn[j]] <- delta.quartet[i]
          break
        }

      }

    }

  }

  #delta.taxa <- as.data.frame(delta.taxa)

  return(list(quartet = delta.quartet, taxa = delta.taxa, taxaNames = taxaNames))

}

#' @title Estimating theta values from a distance matrix.
#'
#' @name esttheta
#'
#' @rdname theta
#'
#' @param disMat a distance matrix with column and row names
#'
#' @return returns a list of two,
#' quartet: vector of theta numbers,
#' each corresponded to a quartet in the distance matrix.
#' taxa: matrix of theta numbers,
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
#' thetas <- esttheta(dis.mat)
#' hist(thetas$quartet)
#'
#' @references
#'
#' @export
esttheta = function(disMat = NULL) {

  taxaNames <- row.names(disMat)
  taxaNumber <- length(taxaNames)

  allCombn <- combn(taxaNumber, 4)

  # for all the combinations  of quartet c(n,4)
  theta.quartet <- vector(mode = "numeric", length = length(allCombn[1,]))


  theta.taxa <- matrix(-1, nr = choose(taxaNumber-1,3), nc = taxaNumber)
  #names(theta.taxa) <- taxaNames

  #theta.quartet <- numeric(length=length(allCombn[1,]))
  for (i in 1:length(allCombn[1,])) {

    aCombn <- allCombn[,i]

    aQuartet <- disMat[aCombn,aCombn]
    aQuartet.melted <- melt(aQuartet)

    aQuartet.melted <- aQuartet.melted[aQuartet.melted$value!=0,]
    sixValues <- aQuartet.melted$value

    threeQuantities <- sort(c(sixValues[3]+sixValues[4],
                              sixValues[2]+sixValues[5],
                              sixValues[1]+sixValues[6]))

    theta.quartet[i] <- (threeQuantities[3] - threeQuantities[2]) / (threeQuantities[3] - threeQuantities[1])

    for (j in 1:4) {

      for (k in 1:length(theta.taxa[,aCombn[j]])) {

        if (theta.taxa[k, aCombn[j]] == -1) {
          theta.taxa[k, aCombn[j]] <- theta.quartet[i]
          break
        }

      }

    }

  }

  #theta.taxa <- as.data.frame(theta.taxa)

  return(list(quartet = theta.quartet, taxa = theta.taxa, taxaNames = taxaNames))

}

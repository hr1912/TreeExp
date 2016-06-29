#' @title Estimating theta values from a distance matrix.
#'
#' @name esttheta
#'
#' @rdname theta
#'
#' @param disMat a distance matrix with column and row names
#'
#' @return returns a vector of theta numbers,
#' each corresponded to a quartet in the distance matrix.
#'
#' @details
#'
#' @export
#'
#' @examples
#' data(tetraexp)
#' dis.mat <- expdist(tetraexp.objects, taxa = "all",
#'                      subtaxa = "Brain",
#'                      method = "rho")
#' thetas <- esttheta(dis.mat)
#' hist(thetas)
#'
#' @references
#'
esttheta = function(disMat = NULL) {

  taxaNames <- row.names(disMat)
  taxaNumber <- length(taxaNames)

  allCombn <- combn(taxaNumber, 4)

  # for all the combinations  of quartet c(n,4)
  theta <- vector(mode = "numeric", length = length(allCombn[1,]))

  for (i in 1:length(allCombn[1,])) {

    aCombn <- allCombn[,i]

    aQuartet <- dis.mat[aCombn,aCombn]
    aQuartet.melted <- melt(aQuartet)

    aQuartet.melted <- aQuartet.melted[aQuartet.melted$value!=0,]
    sixValues <- aQuartet.melted$value

    threeQuantities <- sort(c(sixValues[3]+sixValues[4],
                              sixValues[2]+sixValues[5],
                              sixValues[1]+sixValues[6]))

    theta[i] <- (threeQuantities[3] - threeQuantities[2]) / (threeQuantities[3] - threeQuantities[1])

  }

  theta
}


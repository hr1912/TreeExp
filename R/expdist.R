#‘
#' @title Generate a distance matrix from a taxaExp class
#'
#' @name expdist
#'
#' @param objects a vector of objects of class \code{taxonExp} or an object of class \code{taxaExp}
#' @param taxa one single character or a vector of characters specifying main taxa selected for
#' calculating expression distance.
#' If one single character "all" is given,
#' all the taxa included in the \code{taxaExp} will be matched and selected ("all" by default).
#' @param subtaxa one single character or a vector of characters sepcifying sub taxa selected for
#' calculating expression distance.
#' If one singke character "all" is given,
#' all the subtaxa included in the \code{taxaExp} will be matched and selected ("all" by default).
#' @param rowindex a vector of numbers corresponded to indices of selecting rows
#' @param method specifying which distance method to be used
#' to estimate expression phylogeny in bootstrapping.
#'
#' @return returns an expression distance matrix
#'
#' @examples
#' data(tetraexp)
#' dismat <- expdist(tetraexp.objects, taxa = "all",
#'                  subtaxa = "Brain",
#'                  method = "sou")
#' tr <- root(nj(dismat), "Chicken_Brain")
#' plot(tr)
#'
#' @references
#'
#' @export
expdist = function (objects = NULL, taxa = "all", subtaxa = "all", rowindex = NULL,
                    method = c("sou", "ced", "pea", "souln", "nbdln", "euc", "cos", "jsd"))
{
  #if(verbose) message(date())

  if (is.null(objects) || class(objects) != "taxaExp") {
    stop(paste0(date(), "no valid taxaExp objects input"))
  }

  flag1 <- TRUE
  flag2 <- TRUE

  if (any(grepl("all",taxa, ignore.case = T))) {flag1 = FALSE}
  else { taxa <- gsub("\\s+","",taxa)}

  if (any(grepl("all",subtaxa, ignore.case = T))) {flag2 = FALSE}
  else { subtaxa <- gsub("\\s+","",subtaxa)}

  #browser()
  # subsetting
  if ( flag1 || flag2)

  {

    objects_n <- length(objects)

    objects_new_n <- 0

    #browser()

    for (i in 1:objects_n)

    {
      if (flag1 && flag2) {
        if (any(grepl(objects[[i]]$taxon.name,taxa, ignore.case=T))
            && any(grepl(objects[[i]]$subTaxon.name, subtaxa, ignore.case=T)))
        {objects_new_n <- objects_new_n + 1}

      } else {
        if (any(grepl(objects[[i]]$taxon.name,taxa,ignore.case=T))
            ||  any(grepl(objects[[i]]$subTaxon.name, subtaxa, ignore.case=T)))
        {objects_new_n <- objects_new_n + 1}
      }

    }

    objects_new <- vector("list",length = objects_new_n)

    counter <- 1

    for (i in 1:objects_n)

    {

      if (flag1 && flag2) {
        if (any(grepl(objects[[i]]$taxon.name,taxa,ignore.case=T))
            &&  any(grepl(objects[[i]]$subTaxon.name, subtaxa, ignore.case=T)))
        {
          objects_new[[counter]] <- objects[[i]]
          counter <- counter + 1
        }

      } else {
        if (any(grepl(objects[[i]]$taxon.name,taxa,ignore.case=T))
            ||  any(grepl(objects[[i]]$subTaxon.name, subtaxa, ignore.case=T)))
        {
          objects_new[[counter]] <- objects[[i]]
          counter <- counter + 1
        }
      }

    }

    class(objects_new) <- "taxaExp"

    objects <- objects_new

  }

  #browser()

  method<-match.arg(method)

  message(paste0(date(), ": using ", method, " to calculate pair-wise distance"))

  object_n <- length(objects)

  gene_n <- objects[[1]]$gene.num

  message(paste0(date(),": input ",object_n, " taxa"))
  message(paste0(date(),": total ", gene_n, " genes"))

  #initialization

  #dis.mat <- matrix(0, nr = object_n, nc = object_n)

  reads.count <- matrix(0, nr = gene_n, nc = object_n)
  gene_length <- matrix(0, nr = gene_n, nc = object_n)

  meanRPKM <- matrix(0, nr = gene_n, nc = object_n)

  omega <- vector("numeric", length = object_n)
  taxon.names <- vector("character", length = object_n)

  for (i in 1:object_n) {

    taxon.names[i] = paste0(objects[[i]]$taxon.name, "_", objects[[i]]$subTaxon.name)

    if (is.null(objects[[i]]$omega))
      omega[i] = 0
    else
      omega[i] = objects[[i]]$omega

    gene_length[,i] = objects[[i]]$gene.lengths

    if (is.null(objects[[i]]$readsCount.rmOut))
      reads.count[,i] = apply(objects[[i]]$readsCount.raw,1,mean)
    else
      reads.count[,i] = apply(objects[[i]]$readsCount.rmOut,1,mean)


    if (is.null(objects[[i]]$rpkm.rmOut))
      meanRPKM[,i] = apply(objects[[i]]$rpkm.raw,1,mean)
    else
      meanRPKM[,i] = apply(objects[[i]]$rpkm.rmOut,1,mean)

  }

  if (!is.null(rowindex)) {

    meanRPKM <- meanRPKM[rowindex,]

    reads.count <- reads.count[rowindex,]
    gene_length <- gene_length[rowindex,]

  }

  #browser()
  if (method == "sou")
    dis.mat <- dist.sou(meanRPKM)

  if (method == "ced")
    dis.mat <- dist.ced(meanRPKM)

  if (method == "nbdln")
    dis.mat <- dist.ndbln(reads.count,gene_length,omega)

  if (method == "souln")
    dis.mat <- dist.brownian(reads.count,gene_length)

  if (method == "pea")
    dis.mat <- dist.pea(meanRPKM)

  if (method == "euc")
    dis.mat <- dist.euc(meanRPKM)

  if (method == "cos")
    dis.mat <- dist.cos(meanRPKM)

  if (method == "jsd")
    dis.mat <- dist.jsd(meanRPKM)

  row.names(dis.mat) = taxon.names
  colnames(dis.mat) = taxon.names

  dis.mat

}

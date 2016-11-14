#'
#' @title Expression distance matrix generated from a taxaExp class
#'
#' @name expdist
#' @description Generate An expression distance matrix from a taxaExp class
#' using a specified distance method
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
#' @param logrithm a logical specifying whether to apply expression value log2 tranformation (TRUE by default).
#'
#' @return returns an expression distance matrix
#'
#' @examples
#' data(tetraexp)
#' dismat <- expdist(tetraexp.objects, taxa = "all",
#'                  subtaxa = "Brain",
#'                  method = "pea")
#' tr <- root(NJ(dismat), "Chicken_Brain")
#' plot(tr)
#'
#' @references
#'
#' @export
expdist = function (objects = NULL, taxa = "all", subtaxa = "all", rowindex = NULL,
                    method = c( "pea", "spe","euc", "cos", "jsd",
                                "tani", "jac" ,"u", "nbdln" ), logrithm = TRUE)
{
  #if(verbose) message(date())

  if (is.null(objects) || class(objects) != "taxaExp") {
    stop(paste0(date(), ": no valid taxaExp objects input!"))
  }

  flag1 <- TRUE
  flag2 <- TRUE

  if (any(grepl("all",taxa, ignore.case = T))) {flag1 = FALSE}
  else { taxa <- gsub("\\s+","",taxa)}

  if (any(grepl("all",subtaxa, ignore.case = T))) {flag2 = FALSE}
  else { subtaxa <- gsub("\\s+","",subtaxa)}

  #browser()
  # subsetting
  objects_n <- length(objects)
  objects_new_n <- 0

  if ( flag1 || flag2)

  {
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

  } else {

    objects_new <- vector("list", length = objects_new_n)
    counter <- 1

    for (i in 1:objects_n) {
        objects_new[[counter]] <- objects[[i]]
        counter <- counter + 1
    }

  }

  if (length(objects_new) == 0) {

    stop(paste0(date(),": taxa and subtaxa name not found."))

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

  read.counts <- matrix(0, nr = gene_n, nc = object_n)
  gene_length <- matrix(0, nr = gene_n, nc = object_n)

  expVal <- matrix(0, nr = gene_n, nc = object_n)

  omega <- vector("numeric", length = object_n)
  taxon.names <- vector("character", length = object_n)

  for (i in 1:object_n) {

    taxon.names[i] = paste0(objects[[i]]$taxon.name, "_", objects[[i]]$subTaxon.name)

    if (is.null(objects[[i]]$omega))
      omega[i] = 0
    else
      omega[i] = objects[[i]]$omega

    gene_length[,i] = objects[[i]]$gene.lengths

    if (is.null(objects[[i]]$readCounts.rmOut))
      read.counts[,i] = apply(objects[[i]]$readCounts.raw,1,median)
    else
      read.counts[,i] = apply(objects[[i]]$readCounts.rmOut,1,median)

    expVal[,i] = apply(objects[[i]]$normExp.val,1,median)

  }

  if (!is.null(rowindex)) {

    expVal <- expVal[rowindex,]

    read.counts <- read.counts[rowindex,]
    gene_length <- gene_length[rowindex,]

  }

  if (logrithm) {

    expVal <- apply(expVal, c(1,2), function (x) log2(x+1))

  }
  #browser()

  dis.mat <- switch (method,

    pea = {dist.pea(expVal)},

    spe = {dist.spe(expVal)},

    euc = {dist.euc(expVal)},

    cos = {dist.cos(expVal)},

    jsd = {dist.jsd(expVal)},

    tani = {dist.tani(expVal)},

    jac = {dist.jac(expVal)},

    ced = {dist.ced(expVal)},

    nbdln = {dist.nbdln(read.counts, gene_length, omega)},

    u = {dist.u(read.counts, gene_length)}
  )

  row.names(dis.mat) = taxon.names
  colnames(dis.mat) = taxon.names

  dis.mat + t(dis.mat)
  #as.dist(dis.mat)

}

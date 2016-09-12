

cov.mat = function (meanRPKM = NULL) {

  object_n <- ncol(meanRPKM)
  #gene_n <- nrow(meanRPKM)

  cov.mat <- matrix(0, nr = object_n, nc = object_n)


  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      #V11 <- var(log2(meanRPKM[,i]+1))
      #V22 <- var(log2(meanRPKM[,j]+1))
      #V12 <- cov(log2(meanRPKM[,i]+1), log2(meanRPKM[,j]+1))

      cov.mat[j,i] <- cov(log2(meanRPKM[,i]+1), log2(meanRPKM[,j]+1))

    }

  }

  cov.mat

}



var.arr = function (meanRPKM = NULL) {

  object_n <- ncol(meanRPKM)

  var.array <- vector(mode="numeric", length=object_n)

  for (i in 1:object_n) {

    var.array[i] <- var(log2(meanRPKM[,i]+1))

  }

  var.array

}


estimate.variance = function (objects = NULL, taxa = "all", subtaxa = "all")
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

  #browser()
  covariance.matrix <- cov.mat(meanRPKM)
  variance.array <- var.arr(meanRPKM)

  row.names(covariance.matrix) = taxon.names
  colnames(covariance.matrix) = taxon.names

  names(variance.array) = taxon.names

  return(list(covariance.maxtrix= covariance.matrix, variance.array = variance.array))

}



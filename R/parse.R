#'
#' @title Expression level table from a taxaExp class
#'
#' @name exptabTE
#'
#' @description Generate an expression level table from a taxaExp class
#'
#' @param objects a vector of objects of class \code{taxonExp} or an object of class \code{taxaExp}
#' @param taxa one single character or a vector of characters specifying main taxa to be included in
#' the expression level table.
#' If one single character "all" is given,
#' all the taxa included in the \code{taxaExp} will be matched and included ("all" by default).
#' @param subtaxa one single character or a vector of characters specifying sub taxa to be included in
#' the expression level table.
#' If one single character "all" is given,
#' all the subtaxa included in the \code{taxaExp} will be matched and included ("all" by default).
#' @param rowindex a vector of numbers corresponded to indices of selecting rows
#' @param logrithm a logical specifying whether to apply expression value log2 tranformation (TRUE by default).

#' @return an expression level table: column corresponds to median expression value of all biological samples
#' within one taxa_subtaxa group; row corresponds to othologous genes
#'
#' @export
exptabTE = function (objects = NULL, taxa = "all", subtaxa = "all",
                     rowindex = NULL, logrithm = TRUE)

{

  if (is.null(objects) || class(objects) != "taxaExp") {
    stop(paste0(date(), ": no valid taxaExp object input!"))
  }

  flag1 <- TRUE
  flag2 <- TRUE

  if (any(grepl("all",taxa, ignore.case = T))) {flag1 = FALSE}
  else { taxa <- gsub("\\s+","",taxa)}

  if (any(grepl("all",subtaxa, ignore.case = T))) {flag2 = FALSE}
  else { subtaxa <- gsub("\\s+","",subtaxa)}

  #browser()

  expval_table <- NULL
  sample_names <- NULL
  # subsetting
  objects_n <- length(objects)

  if ( flag1 || flag2)

  {

    for (i in 1:objects_n)

    {
      if (flag1 && flag2) {
        if (any(grepl(objects[[i]]$taxon.name,taxa, ignore.case=T))
            && any(grepl(objects[[i]]$subTaxon.name, subtaxa, ignore.case=T)))
        {
          expval_table <- cbind(expval_table, apply(objects[[i]]$normExp.val, 1, median))
          sample_names <- c(sample_names,
                            paste0(objects[[i]]$taxon.name,"_",objects[[i]]$subTaxon.name))
        }

      } else {
        if (any(grepl(objects[[i]]$taxon.name,taxa,ignore.case=T))
            ||  any(grepl(objects[[i]]$subTaxon.name, subtaxa, ignore.case=T)))
        {
          expval_table <- cbind(expval_table, apply(objects[[i]]$normExp.val, 1, median))
          sample_names <- c(sample_names,
                            paste0(objects[[i]]$taxon.name,"_",objects[[i]]$subTaxon.name))
        }
      }
    }
  } else {

    for (i in 1:objects_n) {
      expval_table <- cbind(expval_table, apply(objects[[i]]$normExp.val, 1, median))
      sample_names <- c(sample_names,
                        paste0(objects[[i]]$taxon.name,"_",objects[[i]]$subTaxon.name))
    }

  }

  if (is.null(expval_table)) {

    stop(paste0(date(),": taxa and subtaxa name not found."))

  } else {

    row.names(expval_table) = objects[[1]]$gene.names
    colnames(expval_table) = sample_names

  }

  if (!is.null(rowindex)) {

    expval_table <- expval_table[rowindex,]

  }

  if (logrithm) {

    expval_table <- apply(expval_table, c(1,2), function(x) log2(x+1))

  }

  expval_table
}


#' @title Expression variance-covariance matrix
#'
#' @name estVarCov
#' @description Generate a variance-covariance matrix from an expression distance matrix
#'
#' @param expMat an expression distance matrix
#'
#' @return returns an expression variance-covariance matrix
#'
#' @export
estVarCov = function (expMat = NULL) {

    object_n <- ncol(expMat)

    cov.mat <- matrix(0, nr = object_n, nc = object_n)

    for (i in 1:(object_n-1)) {

        for (j in (i+1):object_n) {

            cov.mat[j,i] <- cov(expMat[,i], expMat[,j])

        }

    }

    colnames(cov.mat) <- colnames(expMat)
    row.names(cov.mat) <- colnames(expMat)

    cov.mat <- cov.mat + t(cov.mat)

    diag(cov.mat) <- apply(expMat, 2, var)

    cov.mat
}



#' @title Consice display of a taxaExp object
#'
#' @param objects an object of class \code{taxaExp}.
#' @param details a logical specifying whether to print taxa and subtaxa names.
#' @param ... further arguments passed to or from other methods.
#'
#' @return NULL.
#'
#' @examples
#' data(tetraexp)
#' print(tetraexp.objects, details = TRUE)
#'
#' @export
print.taxaExp <- function(objects, details = FALSE, ...) {

  N <- length(objects)
  cat("\n",N, "taxonExp objects", "\n")

  if (details) {
    cat("\n")
      for (i in 1:N) {
        cat ("object", i, ":", objects[[i]]$taxon.name,
             "\t", objects[[i]]$subTaxon.name, "\n")
      }
  }
}

#' Concise display of a taxonExp object
#'
#' @param object an object of class \code{taxonExp}.
#' @param printlen the number of biological replicates title to print (6 by default).
#' @param ... further arguments passed to or from other methods.
#'
#' @return NULL.
#' @export
#'
#' @examples
#' data(tetraexp)
#' print(tetraexp.objects[[1]], printlen = 6)
#'
print.taxonExp <- function(object, printlen = 6, ...) {

  cat("\nOne taxonExp object\n")

  cat("Taxon name: ", object$taxon.name, "\n")

  cat("Subtaxon name: ", object$subTaxon.name, "\n")

  cat("Total gene number: ", object$gene.num, "\n")

  cat("Total bio replicates number: ", object$bioRep.num, "\n")

  cat("Bio replicates titles:\n")
  if (object$bioRep.num > printlen) {
    cat(paste0("\t", paste(object$bioRep.id[1:printlen]),
                collapse = ", "), ", ...\n")

  } else {
    #cat(paste0("\t\t", paste(object$bioRep.id[1:object$bioRep.num]),
    #          collapes = ", "), "\n")
    print(unlist(object$bioRep.id))
  }

  if (is.null(object$readCounts.rmOut)) {
    cat("Outliers NOT removed\n")
  } else {
    cat("Outliers removed\n")
  }

  if (is.null(object$normExp.val)) {
    cat("Normalized expression value NOT calculated\n")
  } else {
    cat("Normalized expression value calculated\n")
    cat("Normalized method: ", object$normalize, "\n")
  }

  if (is.null(object$omega)) {
    cat("Over-dispersion parameter omega NOT calculated\n")
  } else{
    cat("Over-dispersion parameter omega: ", object$omega, "\n")
  }

  cat("\n")
}


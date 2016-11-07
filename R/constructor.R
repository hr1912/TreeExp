
#' @title Construct a taxaExp object.
#'
#' @description  constructor function for \code{taxaExp} objects
#' This function takes in a read counts file and a gene information file.
#' And construct a \code{taxaExp} object from which user can extract information
#' for display or for further analysis.
#'
#' @name TEconstruct
#'
#' @rdname TEconstruct
#'
#' @param readCountsFP a text file contains raw read counts data.
#' Row names correspond with gene names,
#' and column names correspond with taxon and subtaxon names.
#' @param geneInfoFP a text file contains gene length information.
#' Row names conrrespond with gene names,
#' and column names correspond with taxon names.
#' @param taxa one single string or a vector of strings specifying main taxa selected for
#' constructing \code{taxaExp} object.
#' Taxa names are extracted from row names given in gene length file.
#' If one single string "all" is given,
#' all the taxa in the row names will be matched and selected ("all" by default).
#' @param subtaxa one single string or a vector of strings sepcifying sub taxa selected for
#' constructing \code{taxaExp} object.
#' If one single string "all" is given,
#' all the sub taxa in the row names will be matched and selected ("all" by default).
#' @param normalize a single string sepcifying normalization method ("TPM", "RPKM" or "CPM")
#' while constructing \code{taxaExp} object ("TPM" by default).
#' @param rmOut a logical sepcifying whether to remove read counts outliers
#' while constructing \code{taxaExp} objects (TRUE by default).
#' @param verbose a logical specifying whether to print more information on the screen
#' while constructing \code{taxaExp} objects (FALSE by default).
#'
#' @return returns an object of class \code{Taxa} (S3 class, a list of \code{taxonExp} objects).
#'
#' @examples
#'
#' taxa.objects = TEconstruct(readCountsFP = system.file('extdata/tetraexp.read.counts.raw.txt', package = 'TreeExp'),
#'    geneInfoFP = system.file('extdata/tetraexp.length.ortholog.txt', package = 'TreeExp'),
#'    taxa = "all", subtaxa = c("Brain", "Cerebellum"), normalize = "TPM")
#'
#' @export
TEconstruct = function(readCountsFP=NULL, geneInfoFP=NULL, taxa="all", subtaxa="all",
                       normalize=c("TPM", "RPKM", "CPM"), rmOut=FALSE, verbose=FALSE) {

  # check file handle
  if(all(c(is.null(readCountsFP), is.null(geneInfoFP)))){
    stop(paste0(date(),": must provide read counts file path and gene length file path"))
  }

  # check file existance
  if(all(c(!file.exists(readCountsFP), !file.exists(geneInfoFP)))){
    stop(paste0(date(),": fail to open file, check your filename or path"))
  }
  #browser()

  # input
  read.counts.df <- read.table(readCountsFP,header=T)
  row.names(read.counts.df) <- read.counts.df[,1]
  read.counts.df <- read.counts.df[,-1]

  gene.info.df <- read.table(geneInfoFP,header=T)


  if (nrow(read.counts.df) != nrow(gene.info.df)) {
    stop(paste0(date(),": different row length in read countss and gene info files"))
  }

  # remove sample with low read counts

  invalid_arr <- NULL

  for (i in 2:ncol(read.counts.df)) {
    if (mean(read.counts.df[,i]) < 1) {
      invalid_arr <- c(invalid_arr,i)
    }
  }

  message(paste0(date(),": removing ", length(invalid_arr), " sample(s) with ultra-low read counts"))

  if (length(invalid_arr) != 0) {
    invalid_arr = 0 - invalid_arr
    read.counts.df <- read.counts.df[,invalid_arr]
  }

  # gene number and taxon number
  gene_n <- nrow(gene.info.df)

  # get taxon names from read counts file
  taxon_names <- unique(lapply(colnames(read.counts.df), function(x) unlist(strsplit(x, "_"))[1]))
  taxon_n <- length(taxon_names)

  if (taxon_n > ncol(gene.info.df)) {
    stop(paste0(date(),": missing taxa gene length info in gene info file"))
  }

  if (taxon_n < ncol(gene.info.df)) {
    message(paste0(date(),": skipping some taxa gene length info"))
  }

  normalize<-match.arg(normalize)
  message(paste0(date(), ": using ", normalize, " to normalize raw read counts"))
  
  # get taxon names
  #browser()
  cat("\n")
  message(paste0(date(),": start constructiong TE objects"))


  if (!any(grepl("all", taxa, ignore.case = T))) {

    taxon_names <- gsub("\\s+", "", taxa)
    taxon_n <- length(taxon_names)

  }

  message(paste0(date(),": total Taxon number ", taxon_n))

  #browser()

  # get subtaxon number
  subtaxon_names <- unique(lapply(colnames(read.counts.df), function(x) unlist(strsplit(x, "_"))[2]))
  subtaxon_n <- length(subtaxon_names)


  if (!any(grepl("all", subtaxa, ignore.case = T))) {

    subtaxon_names <- gsub("\\s+", "", subtaxa)
    subtaxon_n <- length(subtaxon_names)

  }

  message(paste0(date(),": total sub taxon number ", subtaxon_n))

  title <- lapply(colnames(read.counts.df), function(x) unlist(strsplit(x, "_"))[1]) # taxon names
  subtitle <- lapply(colnames(read.counts.df), function(x) unlist(strsplit(x,"_"))[2]) # subtaxon names
  first_two_names <- unique(paste(title,subtitle,sep="_"))

  index <- intersect(unlist(lapply(taxon_names, function(x) grep(x, first_two_names, ignore.case = T))),
          unlist(lapply(subtaxon_names,function(x) grep(x, first_two_names, ignore.case = T))))

  objects_names <- first_two_names[index]

  objects_number <- length(objects_names)

  #browser()

  cat("\n")
  # get gene names
  #gene.names <- read.counts.df[,1]

  message(paste0(date(),": now constructing ",objects_number, " TE objects..."))

  if (!verbose) progbar <- txtProgressBar(style = 3)

  # initialization

  taxonExp.objects <- vector("list",length = objects_number)
  # the number of TE objects constructed is based on seleted numnber

  # for each taxon

  objects_counter <- 0

  for (i in 1:objects_number) {

    #browser()
    if (verbose) message(paste0(date(),": proceeding taxon ", objects_names[i]))

    # get all the sample names matching objects names
    # bundle all the biological replicates into one TE object

    #browser()

    ttl <- unlist(strsplit(objects_names[i], "_"))[1] #taxon title
    subttl <- unlist(strsplit(objects_names[i], "_"))[2] # subtaxon title
    #ttl <- lapply(names, function(x) unlist(strsplit(x, "_"))[1]) # taxon names
    #subttl <- lapply(names, function(x) unlist(strsplit(x,"_"))[2]) # subtaxon names

    idx <- grep(objects_names[i],colnames(read.counts.df), ignore.case = T)
    names <- strsplit(colnames(read.counts.df)[idx],"_")
    repttl <- unlist(lapply(names, function(x) unlist(strsplit(x,"_"))[3])) # biological replicates title names


    # get gene names and lengths

    gene_info = gene.info.df[grep(ttl,colnames(gene.info.df), ignore.case = T)]
    tmp = apply(gene_info,1,function(x) unlist(strsplit(as.character(x),":")))

    #
    gene_names <- tmp[1,]  # gene names
    gene_lengths <- as.integer(tmp[2,]) # gene lengths


    # foreach subtaxon
    bio_rep_n <- length(repttl) # biological replicates number
    omega <- NULL # omega estimated overdispersion parameter

    read.counts.raw <- apply(read.counts.df[idx], c(1,2), as.numeric)

    read.counts.rmOut <- NULL
    norm_exp_val <- NULL


    if (rmOut) {  # removing outliers loop

      if (verbose) message(paste0(date(),": removing outliers"))

      # remove outliers
      # if read counts is within the quantile(2.5%, 97.5%),
      # remain the original count
      # if read counts is vithin the quantile(97.5%, 100%),
      # replace the original count by the medium count within this quantile

      read.counts.rmOut <- read.counts.raw

      quants.mat<-matrix(data=0,nrow=6,ncol=5) # 0.975 0.5 0.025 0.9875 0.0125

      for (l in 1: bio_rep_n) {
        quants.mat[l,] <- quantile(read.counts.raw[,l],
                                   prob=c(.975,.5,.025,.9875,.0125))
      }

      for (k in 1:gene_n) {

        for (l in 1:bio_rep_n) {

          if (as.numeric(read.counts.rmOut[k,l])>quants.mat[l,1]) {
            read.counts.rmOut[k,l]<-quants.mat[l,4]
          }

          if  (as.numeric(read.counts.rmOut[k,l])<quants.mat[l,3]) {
            read.counts.rmOut[k,l]<-quants.mat[l,5]
          }

        }

      }

    } # removing outliers loop end

    # calculating normalized expression values loop

    if (verbose) message(paste0(date(), ": normalizing raw read counts using", normalize))

    norm_exp_val <- read.counts.raw

    #browser()

    for (l in 1:bio_rep_n) {

      norm_exp_val[,l] <- switch(normalize,

        TPM = { # Transcripts Per kilobase Million
          if (rmOut) {
            read.counts.rmOut[,l] / (gene_lengths / 10^3) /
            (sum(read.counts.rmOut[,l] / (gene_lengths / 10^3)) / 10^6)
          } else {
            read.counts.raw[,l] / (gene_lengths / 10^3) /
            (sum(read.counts.raw[,l] / (gene_lengths / 10^3)) / 10^6)
          }
        },

        RPKM = { # Reads Per Kilobase Million
          if (rmOut) {
            read.counts.rmOut[,l] / (sum(read.counts.rmOut[,l]) / 10^6) /
            (gene_lengths / 10^3)
          } else {
            read.counts.raw[,l] / (sum(read.counts.raw[,l]) / 10^6) /
            (gene_lengths / 10^3)
          }

        },

        CPM = { # Counts Per Million not normalized by gene length
          if (rmOut) {
            read.counts.rmOut[,l] / (sum(read.counts.rmOut[,l]) / 10^6)
          } else {
            read.counts.raw[,l] / (sum(read.counts.raw[,l]) /  10^6)
          }

        }
      )

    }

    objects_counter = objects_counter + 1

    if (verbose) message(paste0(date(),": wrapping up into objects"))

    #browser()
    oneObject <- list(readCounts.raw=read.counts.raw, readCounts.rmOut=read.counts.rmOut,
                      taxon.name = ttl,subTaxon.name = subttl,
                      gene.num = gene_n, gene.names = gene_names, gene.lengths = gene_lengths,
                      bioRep.num = bio_rep_n, bioRep.id = repttl, omega = omega,
                      normalize = normalize, normExp.val = norm_exp_val)

    class(oneObject) <- "taxonExp"

    taxonExp.objects[[objects_counter]] <- oneObject

    #browser()

    if (verbose) message(paste0(date(), ": ", objects_counter, " TE objects constructed"))

    if (verbose) cat("\n")

    if (!verbose) setTxtProgressBar(progbar, objects_counter/objects_number)


  }

  class(taxonExp.objects) <- "taxaExp"

  attr(taxonExp.objects, "taxa") <- unlist(taxon_names)
  attr(taxonExp.objects, "subtaxa") <- unlist(subtaxon_names)

  cat("\n")
  message(date(),": construction complete.")

  taxonExp.objects

}

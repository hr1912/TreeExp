
#' @title Construct a taxaExp object.
#'
#' @description  constructor function for \code{taxaExp} objects
#' This function takes in a reads count file and a gene information file.
#' And construct a \code{taxaExp} object from which user can extract information
#' for display or for further analysis.
#'
#' @name TEconstruct
#'
#' @rdname TEconstruct
#'
#' @param readsCountFP a text file contains raw reads count data.
#' Row names correspond with gene names,
#' and column names correspond with taxon and subtaxon names.
#' @param geneInfoFP a text file contains gene length information.
#' Row names conrrespond with gene names,
#' and column names correspond with taxon names.
#' @param taxa one single character or a vector of characters specifying main taxa selected for
#' constructing \code{taxaExp} object.
#' Taxa names are extracted from row names given in gene length file.
#' If one single character "all" is given,
#' all the taxa in the row names will be matched and selected ("all" by default).
#' @param subtaxa one single character or a vector of characters sepcifying sub taxa selected for
#' constructing \code{taxaExp} object.
#' If one single character "all" is given,
#' all the sub taxa in the row names will be matched and selected ("all" by default).
#' @param calRPKM a logical sepcifying whether to calculate RPKM value
#' while constructing \code{taxaExp} object (TRUE by default).
#' @param rmOut a logical sepcifying whether to remove reads count outliers
#' while constructing \code{taxaExp} objects (TRUE by default).
#' @param verbose a logical specifying whether to print more information on the screen
#' while constructing \code{taxaExp} objects (FALSE by default).
#'
#' @return returns an object of class \code{Taxa} (S3 class, a list of \code{taxonExp} objects).
#'
#' @export
#'
#' @examples
#'
#' taxa.objects = TEconstruct(readsCountFP = system.file('extdata/tetraexp.reads.count.raw.txt', package = 'phyExp'),
#'    geneInfoFP = system.file('extdata/tetraexp.gene.length.ortholog.txt', package = 'phyExp'),
#'    taxa = "all", subtaxa = c("Brain", "Cerebellum"), calRPKM = TRUE, rmOut =TRUE)
#'
#'
TEconstruct = function(readsCountFP=NULL, geneInfoFP=NULL, taxa="all", subtaxa="all",
                       calRPKM=TRUE, rmOut=TRUE, verbose=FALSE) {

  # check file handle
  if(all(c(is.null(readsCountFP), is.null(geneInfoFP)))){
    stop(paste0(date(),": must provide reads count file path and gene length file path"))
  }

  # check file existance
  if(all(c(!file.exists(readsCountFP), !file.exists(geneInfoFP)))){
    stop(paste0(date(),": fail to open file, check your filename or path"))
  }


  #browser()

  # input
  reads.count.df <- read.table(readsCountFP,header=T)
  row.names(reads.count.df) <- reads.count.df[,1]
  reads.count.df <- reads.count.df[,-1]

  gene.info.df <- read.table(geneInfoFP,header=T)


  if (nrow(reads.count.df) != nrow(gene.info.df)) {
    stop(paste0(date(),": different row length in reads counts and gene info files"))
  }

  # gene number and taxon number
  gene_n <- nrow(gene.info.df)

  # get taxon names from reads count file
  taxon_names <- unique(lapply(colnames(reads.count.df), function(x) unlist(strsplit(x, "_"))[1]))
  taxon_n <- length(taxon_names)

  if (taxon_n > ncol(gene.info.df)) {
    stop(paste0(date(),": missing taxa gene length info in gene info file"))
  }

  if (taxon_n < ncol(gene.info.df)) {
    message(paste0(date(),": skipping some taxa gene length info"))
  }

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
  subtaxon_names <- unique(lapply(colnames(reads.count.df), function(x) unlist(strsplit(x, "_"))[2]))
  subtaxon_n <- length(subtaxon_names)


  if (!any(grepl("all", subtaxa, ignore.case = T))) {

    subtaxon_names <- gsub("\\s+", "", subtaxa)
    subtaxon_n <- length(subtaxon_names)

  }

  message(paste0(date(),": total sub taxon number ", subtaxon_n))

  title <- lapply(colnames(reads.count.df), function(x) unlist(strsplit(x, "_"))[1]) # taxon names
  subtitle <- lapply(colnames(reads.count.df), function(x) unlist(strsplit(x,"_"))[2]) # subtaxon names
  first_two_names <- unique(paste(title,subtitle,sep="_"))

  index <- intersect(unlist(lapply(taxon_names, function(x) grep(x, first_two_names, ignore.case = T))),
          unlist(lapply(subtaxon_names,function(x) grep(x, first_two_names, ignore.case = T))))

  objects_names <- first_two_names[index]

  objects_number <- length(objects_names)

  #browser()

  cat("\n")
  # get gene names
  #gene.names <- reads.count.df[,1]

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

    idx <- grep(objects_names[i],colnames(reads.count.df), ignore.case = T)
    names <- strsplit(colnames(reads.count.df)[idx],"_")
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

    reads.count.raw <- reads.count.df[idx]

    reads.count.rmOut <- NULL
    rpkm.raw <- NULL
    rpkm.rmOut <- NULL


    if (rmOut) {  # removing outliers loop

      if (verbose) message(paste0(date(),": removing outliers"))

      # remove outliers
      # if reads count is within the quantile(2.5%, 97.5%),
      # remain the original count
      # if reads count is vithin the quantile(97.5%, 100%),
      # replace the original count by the medium count within this quantile

      reads.count.rmOut <- reads.count.raw

      quants.mat<-matrix(data=0,nrow=6,ncol=5) # 0.975 0.5 0.025 0.9875 0.0125

      for (l in 1: bio_rep_n) {
        quants.mat[l,] <- quantile(reads.count.raw[,l],
                                   prob=c(.975,.5,.025,.9875,.0125))
      }

      for (k in 1:gene_n) {

        for (l in 1:bio_rep_n) {

          if (as.numeric(reads.count.rmOut[k,l])>quants.mat[l,1]) {
            reads.count.rmOut[k,l]<-quants.mat[l,4]
          }

          if  (as.numeric(reads.count.rmOut[k,l])<quants.mat[l,3]) {
            reads.count.rmOut[k,l]<-quants.mat[l,5]
          }

        }

      }

    } # removing outliers loop end

    # calculating RPKM values loop

    if (calRPKM) {

      rpkm.raw <- reads.count.raw

      if (verbose) message(paste0(date(),": calculating RPKM from raw reads count"))

      if (rmOut) {
        if (verbose) message(paste0(date(),": calculating RPKM from outlier removed reads count"))
        rpkm.rmOut <- rpkm.raw
      }

      for (l in 1:bio_rep_n) {

        total_reads_count_raw = sum(reads.count.raw)
        if (rmOut) total_reads_count_rmOut = sum(reads.count.rmOut)

        #browser()

        for (k in 1:gene_n) {
          rpkm.raw[k,l] = reads.count.raw[k,l]/ ( total_reads_count_raw  / 10^6 * gene_lengths[k] / 10^3)
          if (rmOut) rpkm.rmOut[k,l] = reads.count.rmOut[k,l] / (total_reads_count_rmOut / 10^6 * gene_lengths[k] / 10^3)

        }
      }

    } #calculating RPKM values loop end

    objects_counter = objects_counter + 1

    if (verbose) message(paste0(date(),": wrapping up into objects"))

    #browser()
    oneObject <- list(readsCount.raw=reads.count.raw, readsCount.rmOut=reads.count.rmOut,
                      taxon.name = ttl,
                      subTaxon.name = subttl,
                      gene.num = gene_n, gene.names = gene_names, gene.lengths = gene_lengths,
                      bioRep.num = bio_rep_n, bioRep.id = repttl, omega = omega,
                      rpkm.raw = rpkm.raw, rpkm.rmOut = rpkm.rmOut)

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

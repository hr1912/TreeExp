

#' @title Bootstrapping expression phylogeny
#'
#' @description  bootstrap by resampling gene (gene, transcript, exon, etc..)
#'
#' @name boot.exphy
#'
#' @param phy an object of class \code{phylo}.
#' @param objects a vector of objects of class \code{taxonExp} or an object of class \code{taxaExp}.
#' @param rowindex a vector of numbers corresponded to indices of selecting rows
#' @param method  specifying which distance method to be used
#' to estimate expression phylogeny in bootstrapping.
#' @param B the number of bootstrap replicates.
#' @param rooted if "phy" is a rooted tree, a character of the root node's label when constructing "phy";
#' if "phy" is unrooted tree, NULL (NULL by default).
#' @param trees a logical specifying whether to return the bootstrapped trees (FALSE by default).
#'
#' @return similar to boot.phylo in ape, boot.exptree returns a numeric vector
#' which ith element is the number associated to the ith node of phy.
#' if trees = TRUE, boot.exptree returns a list whose first element (named "BP") is like before,
#' and the second element("trees") is a list with the bootstrapped trees.
#'
#' @author Hang Ruan (hang.ruan@hotmail.com).
#'
#' @rdname boot.exphy
#'
#' @examples
#'
#' data(tetraexp)
#' dismat <- expdist(tetraexp.objects, taxa = "all",
#'                  subtaxa = "Brain",
#'                  method = "pea")
#' tr <- root(NJ(dismat), "Chicken_Brain")
#' plot(tr)
#' bs <- boot.exphy(tr, tetraexp.objects, method = "sou",
#'                  B = 100, rooted = "Chicken_Brain")
#' nodelabels(bs)
#'
#' @export
boot.exphy = function (phy = NULL, objects = NULL, rowindex = NULL,
                       method = c( "sou", "pea", "spe","euc", "cos", "jsd",
                                   "tani", "jac" ,"u", "nbdln" ),
                         B = 100, rooted = NULL, trees = FALSE)
{

  method<- match.arg(method)

  message(paste0(date(), ": start bootstapping ", B,  " times using ", method))

  objects_sub_n = length(phy$tip.label)
  object_n = length(objects)

  objects.sub <- vector("list",length = objects_sub_n)

  counter <- 1

  for (i in 1:object_n) {

    if (any(grepl(objects[[i]]$subTaxon.name,phy$tip.label,ignore.case=T))) {
      objects.sub[[counter]] <- objects[[i]]
      counter <- counter + 1
    }

  }

  class(objects.sub) <- "taxaExp"

  #if (counter != (objects_sub_n+1)) {
  #  message(paste0(date(),"incompatible phylo object and taxaExp objects"))
  #  stop()
  #}

  objects_sub_n <- length(objects.sub)

  gene_n <- objects.sub[[1]]$gene.num

  message(paste0(date(),": input ", objects_sub_n, " taxa"))
  message(paste0(date(),": total ", gene_n, " genes"))

  read.counts <- matrix(0, nr = gene_n, nc = objects_sub_n)
  gene_length <- matrix(0, nr = gene_n, nc = objects_sub_n)

  expVal <- matrix(0, nr = gene_n, nc = objects_sub_n)

  rmOut.flag = 1
  if (is.null(objects.sub[[1]]$readCounts.rmOut)) {rmOut.flag = 0}

  taxon.names <- vector("character", length = objects_sub_n)

  for (i in 1:objects_sub_n) {

    taxon.names[i] = paste0(objects.sub[[i]]$taxon.name, "_", objects.sub[[i]]$subTaxon.name)

    if (rmOut.flag) {
      read.counts[,i] = apply(objects.sub[[i]]$readCounts.rmOut,1,median)
    } else {
      read.counts[,i] = apply(objects.sub[[i]]$readCounts.raw,1,median)
    }

    gene_length[,i] = objects.sub[[i]]$gene.lengths

    expVal[,i] = apply(objects.sub[[i]]$normExp.val,1,median)

  }

  if (!is.null(rowindex)) {

    read.counts = read.counts[rowindex,]
    gene_length = gene_length[rowindex,]

    expVal = expVal[rowindex,]

  }

  boot.tree <- vector("list",B)

  message(paste0(date(),": calculating bootstrap values..."))

  progbar <- txtProgressBar(style = 3)

  if (is.null(rowindex)) {
    y <- gene_n
  } else {
    y <- length(rowindex)
  }

  for (n in 1:B) {

    gene_index <- unlist(sample(y, replace = T))

    expVal.samp <- expVal[gene_index,]

    dis.mat <- switch(method,

      sou = {dist.sou(expVal.samp)},

      pea = {dist.pea(expVal.samp)},

      spe = {dist.spe(expVal.samp)},

      euc = {dist.euc(expVal.samp)},

      cos = {dist.cos(expVal.samp)},

      jsd = {dist.jsd(expVal.samp)},

      tani = {dist.tani(expVal.samp)},

      jac = {dist.jac(expVal.samp)},

      u = {

        read.counts.samp <- read.counts[gene_index,]
        gene_length.samp <- gene_length[gene_index,]

        dist.u(read.counts.samp,gene_length.samp)

      },

      nbdln = {

        read.counts.samp <- read.counts[gene_index,]
        gene_length.samp <- gene_length[gene_index,]
        omega.samp <- .estomega.sample(objects.sub,gene_index)

        dist.nbdln(read.counts.samp,gene_length.samp,omega.samp)
      }

    )

    row.names(dis.mat) = taxon.names
    colnames(dis.mat) = taxon.names

    #browser()

    if (!is.null(rooted))
      boot.tree[[n]] <- root(nj(dis.mat),rooted)
    else
      boot.tree[[n]] <- nj(dis.mat)


    setTxtProgressBar(progbar, n/B)

  }

  close(progbar)

  message(paste0(date(),": calculating bootstrap values..."))

  #browser()

  pp <- prop.part(boot.tree)
  ans <- prop.clades(phy, part = pp, rooted = !is.null(rooted))

  if (trees) {
    class(boot.tree) <- "multiPhylo"
    ans <- list(BP = ans, trees = boot.tree)
  }

  message(paste0(date(), ": done bootstrapping"))

  ans

}

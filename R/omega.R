
#' @title Estimation of over-dispersion parameter omega
#'
#' @param objects a vector of objects of class \code{taxonExp} or an object of class \code{taxaExp}
#' @param overwrite a logical specifying whether to overwrite the existing omega value
#' stored in the \code{taxonExp} objects
#'
#' @name estomega
#'
#' @rdname estomega
#'
#' @examples
#' data(tetraexp)
#' tetraexp.objects.with.omega <- estomega(tetraexp.objects)
#' tetraexp.objects.with.omega[[1]]$omega


estomega = function (objects = NULL, overwrite = TRUE) {

  object_n = length(objects)

  message(paste0(date(), ": starting omega estimation"))
  message(paste0(date(),": input ",object_n, " taxa"))

  objects_add_omega = objects

  for (i in 1:object_n) {

    #
    omega = objects[[i]]$omega
    bio_rep_n = objects[[i]]$bioRep.num

    # omega already calculated or no biological replicates

    #browser()

    if ((is.null(omega) || bio_rep_n ==1) && !overwrite ) next

    reads.count.rmOut = objects[[i]]$readsCount.rmOut

    if (is.null(reads.count.rmOut)) {
      message("\n no outlier removed reads count data for object ", i,
          "which is a prerequisite for omega estimation, skipping ...\n")
      next
    }

    gene_n = objects[[i]]$gene.num
    gene_length = objects[[i]]$gene.length

    gene_length_sum = sum(gene_length)

    rk_var<-rep(0,bio_rep_n)

    for (j in 1:bio_rep_n) {
      rk_var[j] <- sum(reads.count.rmOut[,j]) / gene_length_sum
    }

    r0_var <- mean(rk_var)
    rstar_var <- rk_var^2
    rstar_var <- sum(rstar_var) / length(rstar_var)

    # calculating variable a and b for each biological replicates

    a<- rstar_var/rk_var^2
    b<- r0_var/rk_var


    reads_count_exp <-rep(0,gene_n)
    reads_count_square_exp<-rep(0,gene_n)
    reads_count_exp_square<-rep(0,gene_n)
    reads_count_var<-rep(0,gene_n)


    for (k in 1:gene_n) {

      for (j in 1:bio_rep_n) {

        reads_count_exp[k] <- reads_count_exp[k] + b[j] * as.numeric(reads.count.rmOut[k,j])
        reads_count_square_exp[k] <- reads_count_square_exp[k] +  a[j] * (as.numeric(reads.count.rmOut[k,j]))^2

      }

      reads_count_exp[k] = reads_count_exp[k] / bio_rep_n  # m_i mean of reads count for each gene
      reads_count_square_exp[k] = reads_count_square_exp[k] / bio_rep_n  #

      reads_count_var[k] = reads_count_square_exp[k] - reads_count_exp[k]^2 # V_i variance of reads count for each gene E[x2i] - (xi)2
      reads_count_exp_square[k] = reads_count_exp[k]^2


    }

    sum_vi_var<-sum(reads_count_var)
    #sum_mi_var<-sum(reads_count_exp)

    #sum_mi_square_var <- sum(reads_count_square_exp)
    sum_mi_square <- sum(reads_count_exp_square)


    #omega<-(sum_vi_var-sum_mi_var) * gene_n / sum_mi_var^2
    omega<- sum_vi_var / sum_mi_square / bio_rep_n


    objects_add_omega[[i]]$omega = omega

  }

  message(paste0(date(),": omega estimation finished for ", object_n, " objects"))

  return(objects_add_omega)

}

#' @title Internal function for omega estimation in bootstrapping
#'
#' @name estomega.sample
#'
estomega.sample = function (objects = NULL, geneIndex = NULL) {

  objects_n = length(objects)

  omega <- vector("numeric", length = objects_n)

  for (i in 1:objects_n) {

    reads.count.rmOut = objects[[i]]$readsCount.rmOut[geneIndex,]

    bio_rep_n = objects[[i]]$bioRep.num
    gene_n = objects[[i]]$gene.num

    if (bio_rep_n < 2) { omega[i] = 0; next }

    gene_length = objects[[i]]$gene.lengths[geneIndex]

    gene_length_sum = sum(gene_length)
    rk_var<-rep(0,bio_rep_n)

    for (j in 1:bio_rep_n) {
      rk_var[j] <- sum(reads.count.rmOut[,j]) / gene_length_sum
    }

    r0_var <- mean(rk_var)
    rstar_var <- rk_var^2
    rstar_var <- sum(rstar_var) / length(rstar_var)

    # calculating variable a and b for each biological replicates

    a<- rstar_var/rk_var^2
    b<- r0_var/rk_var


    reads_count_exp <-rep(0,gene_n)
    reads_count_square_exp<-rep(0,gene_n)
    reads_count_exp_square<-rep(0,gene_n)
    reads_count_var<-rep(0,gene_n)

    for (k in 1:gene_n) {

      for (j in 1:bio_rep_n) {

        reads_count_exp[k] <- reads_count_exp[k] + b[j] * as.numeric(reads.count.rmOut[k,j])
        reads_count_square_exp[k] <- reads_count_square_exp[k] +  a[j] * (as.numeric(reads.count.rmOut[k,j]))^2

      }

      reads_count_exp[k] = reads_count_exp[k] / bio_rep_n
      reads_count_square_exp[k] = reads_count_square_exp[k] / bio_rep_n

      reads_count_var[k] = reads_count_square_exp[k] - reads_count_exp[k]^2
      reads_count_exp_square[k] = reads_count_exp[k]^2


    }

    #sum_mi_square_var <- sum(reads_count_square_exp)
    sum_mi_square <- sum(reads_count_exp_square)

    sum_vi_var<-sum(reads_count_var)
    #sum_mi_var<-sum(reads_count_exp)

    #omega<-(sum_vi_var-sum_mi_var) * gene_n / sum_mi_var^2
    omega[i]<- sum_vi_var / sum_mi_square / bio_rep_n

  }

  return(omega)

}


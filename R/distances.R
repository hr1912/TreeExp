
# Internal function for estimating Euclidean distance
.euc.dist = function(x,y) sqrt(sum((x-y)^2))

# Internal function for estimating Cosine similarity
.cosine.sim = function(x,y) sum(x*y) / (sqrt(sum(x^2))*sqrt(sum(y^2)))

# Internal function for estimating Kullback-Leibler Divergence
.kld = function(x,y) sum(x * log2(x/y))

# Internal function for estimating Tanimoto distance
.tani.dist = function(x,y) sum(pmax(x,y) - pmin(x,y)) / sum(pmax(x,y))

# Internal function for estimating Jaccard similarity
.jac.sim = function(x,y) sum(x*y) / (sum(x^2) + sum(y^2) - sum(x*y))

#' @title Internal functions for estimating pair-wise expression distances
#'
#' @name distances
#'
#' @rdname distances
#'
#' @description Several published methods to estimate pair-wise expression distances
#'
#' @references
#' Chen H, He X. 2016. The Convergent Cancer Evolution toward a Single Cellular Destination.
#' Mol Biol Evol 33:4-12
#'
#' Gu X, Su Z. 2007. Tissue-driven hypothesis of genomic evolution and sequence-expression correlations.
#' Proc Natl Acad Sci USA 104:2779-2784.
#'
#' Pereira V, Waxman D, Eyre-Walker A. 2009. A problem with the correlation coefficient as a measure of gene expression divergence.
#' Genetics 183:1597-1600.
#'
#' Sudmant PH, Alexis MS, Burge CB. 2015. Meta-analysis of RNA-seq expression data across species, tissues and studies.
#' Genome Biol 16:287.
#'
#' Distance based on negative bio distribution and log normal model
#' @export
dist.nbdln = function (reads.count = NULL, gene_length = NULL, omega = NULL) {

  object_n <- ncol(reads.count)
  gene_n <- nrow(reads.count)

  dis.mat <- matrix(0, nr = object_n, nc = object_n)

  rk_var <- rep(0, object_n)

  ci <- vector("list", length = object_n)

  #Rx <- rep(0, object_n)

  for (i in 1:object_n) {

    gene_length_sum = sum(gene_length[,i])

    rk_var[i] = sum(reads.count[,i]) / gene_length_sum

    ci[[i]] <- gene_length[,i] / mean(gene_length[,i])
    #Rx[i] <- sum(reads.count[,i]) / gene_n
  }

  #R0 <- mean(Rx)

  # calculate variable c

  #browser()

  #gene_length_square <- apply(gene_length,c(1,2), function(x) x^2)

  #sum_gene_length <- apply(gene_length,2, sum)
  #sum_gene_length_square <- apply(gene_length_square,2,sum)

  #c_var <- gene_n * sum_gene_length_square / sum_gene_length ^2


  # calculate variable a and b

  r0_var <- mean(rk_var)
  rstar_var <- rk_var^2
  rstar_var <- sum(rstar_var) / length(rstar_var)

  a_var <- rstar_var / rk_var^2
  b_var <- r0_var / rk_var

  # calculate variable B
  #Bx <- Rx / R0

  #browser()

  for (i in 1:(object_n-1)) {

    #browser()

    JXX = a_var[i]/(sum(ci[[i]]^2)/gene_n)*sum(as.numeric(lapply(reads.count[,i], function(x) x^2 )))/gene_n
    - b_var[i]*mean(reads.count[,i])

    #JXX = a_var[i]/c_var[i]*sum(as.numeric(lapply(reads.count[,i], function(x) x^2 )))/gene_n
    #- b_var[i]*mean(reads.count[,i])

    #JXX = sum(as.numeric(lapply(reads.count[,i], function(x) x^2 )))/gene_n/(Bx[i]^2*c_var[i])
    #  - mean(reads.count[,i])/Bx[i]

    #browser()

    for (j in (i+1):object_n) {

      JYY = a_var[j]/(sum(ci[[j]]^2)/gene_n)*sum(as.numeric(lapply(reads.count[,j], function(x) x^2 )))/gene_n
      - b_var[j]*mean(reads.count[,j])

      #JYY = a_var[j]/c_var[j]*sum(as.numeric(lapply(reads.count[,j], function(x) x^2 )))/gene_n
      #- b_var[j]*mean(reads.count[,i])

      #JYY = sum(as.numeric(lapply(reads.count[,j], function(x) x^2 )))/gene_n/(Bx[j]^2*c_var[j])
      #  - mean(reads.count[,j])/Bx[j]

      JXY = b_var[i] * b_var[j] / (sum(ci[[i]]*ci[[j]])/gene_n) * mean(reads.count[,i]*reads.count[,j])
      #JXY = b_var[i] * b_var[j] / sqrt(c_var[i] * c_var[j]) * mean(reads.count[,i]*reads.count[,j])
      #JXY = mean(reads.count[,i]*reads.count[,j])/(Bx[i]*Bx[j]*c_var[i])

      dis.mat[j,i] = -log(JXY^2 / (JXX*JYY)) - log(1+omega[i]) - log(1+omega[j])
      #dis.mat[j,i] = -log(JXY^2 / (JXX*JYY))
    }
  }


  dis.mat

}

# distance based on stationary OU and Log normal model distance
#' @rdname distances
#'
#' @export
dist.u = function (reads.count = NULL, gene_length = NULL) {


  object_n <- ncol(reads.count)
  gene_n <- nrow(reads.count)

  dis.mat <- matrix(0, nr = object_n, nc = object_n)

  rk_var <- rep(0, object_n)

  ci <- vector("list", length = object_n)

  for (i in 1:object_n) {

    gene_length_sum = sum(gene_length[,i])

    rk_var[i] = sum(reads.count[,i]) / gene_length_sum

    ci[[i]] <- gene_length[,i] / mean(gene_length[,i])
  }

  # calculate variable c

  #browser()

  #gene_length_square <- apply(gene_length,c(1,2), function(x) x^2)

  #sum_gene_length <- apply(gene_length,2, sum)
  #sum_gene_length_square <- apply(gene_length_square,2,sum)

  #c_var <- gene_n * sum_gene_length_square / sum_gene_length ^2

  # calculate variable a and b

  r0_var <- mean(rk_var)
  rstar_var <- rk_var^2
  rstar_var <- sum(rstar_var) / length(rstar_var)

  a_var <- rstar_var / rk_var^2
  b_var <- r0_var / rk_var

  #browser()

  for (i in 1:(object_n-1)) {

    #browser()

    JXX = a_var[i]/(sum(ci[[i]]^2)/gene_n)*sum(as.numeric(lapply(reads.count[,i], function(x) x^2 )))/gene_n
    - b_var[i]*mean(reads.count[,i])
    #JXX = a_var[i]/c_var[i]*sum(as.numeric(lapply(reads.count[,i], function(x) x^2 )))/gene_n
    #- b_var[i]*mean(reads.count[,i])

    #mX = reads.count.mean[i]
    mX = mean(reads.count[,i] * b_var[i] / (sum(ci[[i]]^2/gene_n)))
    #mX = mean(reads.count[,i] * b_var[i] / c_var[i])
    #browser()

    for (j in (i+1):object_n) {

      JYY = a_var[j]/(sum(ci[[j]]^2)/gene_n)*sum(as.numeric(lapply(reads.count[,j], function(x) x^2 )))/gene_n
      - b_var[j]*mean(reads.count[,j])
      #JYY = a_var[j]/c_var[j]*sum(as.numeric(lapply(reads.count[,j], function(x) x^2 )))/gene_n
      #- b_var[j]*mean(reads.count[,j])

      JXY = b_var[i] * b_var[j] / (sum(ci[[i]]*ci[[j]])/gene_n) * mean(reads.count[,i]*reads.count[,j])

      #JXY = b_var[i] * b_var[j] / sqrt(c_var[i] * c_var[j] ) * mean(reads.count[,i]*reads.count[,j])

      #mY = reads.count.mean[j]
      mY = mean(reads.count[,j] * b_var[j] / (sum(ci[[j]]^2/gene_n)))
      #mY = mean(reads.count[,j] * b_var[j] / c_var[j])

      corr = log(JXY/(mX * mY)) / sqrt(log(JXX/mX^2)*log(JYY/mY^2))

      #browser()
      dis.mat[j,i] = -log(corr)
      #dis.mat[j,i] = 1 - corr

    }
  }

  dis.mat

}

# Jesen-Shannon divergence
#' @rdname distances
#'
#' @export
dist.jsd = function (expMat = NULL) {

  object_n <- ncol(expMat)
  gene_n <- nrow(expMat)

  dis.mat <- matrix(0, nr = object_n, nc = object_n)

  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      if (any(expMat[,i] == 0) || any(expMat[,j] == 0)) {

        mu <- (expMat[,i] + expMat[,j] + .02) / 2
        dis.mat[j,i] <- sqrt(.5 * .kld(expMat[,i]+.01, mu) + .5 * .kld(expMat[,j]+.01, mu))

      } else {

        mu <- (expMat[,i] + expMat[,j]) / 2
        dis.mat[j,i] <- sqrt(.5 * .kld(expMat[,i], mu) + .5 * .kld(expMat[,j], mu))

      }

    }

  }

  dis.mat

}

# Pearson distance
#' @rdname distances
#'
#' @export
dist.pea = function (expMat = NULL) {

  object_n <- ncol(expMat)
  gene_n <- nrow(expMat)

  dis.mat <- matrix(0, nr = object_n, nc = object_n)


  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      dis.mat[j,i] <- 1 - cor(expMat[,i],expMat[,j])

    }

  }

  #browser()

  dis.mat

}

# Spearman distance
#' @rdname distances
#'
#' @export
dist.spe = function (expMat = NULL) {

  object_n <- ncol(expMat)
  gene_n <- nrow(expMat)

  dis.mat <- matrix(0, nr = object_n, nc = object_n)


  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      dis.mat[j,i] <- 1 - cor(expMat[,i],expMat[,j],
                              method = "spearman")

    }

  }

  #browser()

  dis.mat

}

# Euclidean distance
#' @rdname distances
#'
#' @export
dist.euc = function (expMat = NULL) {


  object_n <- ncol(expMat)
  gene_n <- nrow(expMat)

  dis.mat <- matrix(0, nr = object_n, nc = object_n)

  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      dis.mat[j,i] <- .euc.dist(expMat[,i],expMat[,j])
    }

  }


  return (dis.mat)

}

# Cosine distance
#' @rdname distances
#'
#' @export dist.cos
dist.cos = function (expMat = NULL) {

  object_n <- ncol(expMat)
  gene_n <- nrow(expMat)

  dis.mat <- matrix(0, nr = object_n, nc = object_n)


  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      dis.mat[j,i] <- 1-.cosine.sim(expMat[,i],expMat[,j])

    }

  }

  dis.mat

}

# Tanimoto distance
#' @rdname distances
#'
#' @export
dist.tani = function (expMat = NULL) {

  object_n <- ncol(expMat)
  gene_n <- nrow(expMat)

  dis.mat <- matrix(0, nr = object_n, nc = object_n)


  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      dis.mat[j,i] <- .tani.dist(expMat[,i],expMat[,j])

    }

  }


  dis.mat

}

# Jaccard distance
#' @rdname distances
#'
#' @export
dist.jac = function (expMat = NULL) {

  object_n <- ncol(expMat)
  gene_n <- nrow(expMat)

  dis.mat <- matrix(0, nr = object_n, nc = object_n)


  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      dis.mat[j,i] <- 1 - .jac.sim(expMat[,i],expMat[,j])

    }

  }


  dis.mat

}

# Converntional expression distance
#' @rdname distances
#'
#' @export
dist.ced = function (expMat = NULL) {

  object_n <- ncol(expMat)
  gene_n <- nrow(expMat)

  dis.mat <- matrix(0, nr = object_n, nc = object_n)


  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      #dis.mat[j,i] <- V11+V22-2*V12
      dis.mat[j,i] <- (.euc.dist(expMat[,i],expMat[,j]))^2 / gene_n

    }

  }

  dis.mat

}

# Distance based on stationary Ornstein-Uhlenback model
#
#' @rdname distances
#'
#' @export
dist.sou = function (expMat = NULL) {

  object_n <- ncol(expMat)
  #gene_n <- nrow(expMat)

  dis.mat <- matrix(0, nr = object_n, nc = object_n)


  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      V11 <- var(expMat[,i])
      V22 <- var(expMat[,j])
      V12 <- cov(expMat[,i], expMat[,j])

      dis.mat[j,i] <- -log(V12/sqrt(V11*V22))

    }

  }

  dis.mat

}


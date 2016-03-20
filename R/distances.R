
# Internal function for estimating Euclidean distance
euc.dist = function(x,y) sqrt(sum((x-y)^2))

# Internal function for estimating Cosine similarity
cosine.sim = function (x,y) sum(x*y) / (sqrt(sum(x^2))*sqrt(sum(y^2)))

# Internal function for estimating Kullback-Leibler Divergence
kld = function (x,y) sum(x * log2(x/y))

#' @title Internal functions for estimating pair-wise expression distances
#'
#' @name distances
#'
#' @rdname distances
#'
#' @description Several published methods to estimate pair-wise expression distances
#'
#' @references
#'
#' Distance based on negative bio distribution and log normal model
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
dist.souln = function (reads.count = NULL, gene_length = NULL) {


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
dist.jsd = function (meanRPKM = NULL, taxon.names = NULL) {

  object_n <- ncol(meanRPKM)
  gene_n <- nrow(meanRPKM)

  dis.mat <- matrix(0, nr = object_n, nc = object_n)

  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      mu <- (meanRPKM[,i]+1 + meanRPKM[,j]+1) / 2
      dis.mat[j,i] <- sqrt(.5 * kld(meanRPKM[,i]+1, mu) + .5 * kld(meanRPKM[,j]+1, mu))

    }

  }

  dis.mat

}

# 1-rho distance
#' @rdname distances
dist.rho = function (meanRPKM = NULL) {

  object_n <- ncol(meanRPKM)
  gene_n <- nrow(meanRPKM)

  dis.mat <- matrix(0, nr = object_n, nc = object_n)


  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      dis.mat[j,i] <- 1 - cor(log2(meanRPKM[,i]+1),log2(meanRPKM[,j]+1),method = "spearman")

    }

  }

  #browser()

  dis.mat


}

# Euclidean distance
#' @rdname distances
dist.euc = function (meanRPKM = NULL) {


  object_n <- ncol(meanRPKM)
  gene_n <- nrow(meanRPKM)

  dis.mat <- matrix(0, nr = object_n, nc = object_n)

  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      #dis.mat[j,i] <- (euc.dist(log2(meanRPKM[,i]+1),log2(meanRPKM[,j]+1)))/gene_n
      dis.mat[j,i] <- (euc.dist(log2(meanRPKM[,i]+1),log2(meanRPKM[,j]+1)))
    }

  }


  return (dis.mat)

}

# Cosine distance
#' @rdname distances
dist.cos = function (meanRPKM = NULL) {

  object_n <- ncol(meanRPKM)
  gene_n <- nrow(meanRPKM)

  dis.mat <- matrix(0, nr = object_n, nc = object_n)


  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      dis.mat[j,i] <- 1-cosine.sim(log2(meanRPKM[,i]+1),log2(meanRPKM[,j]+1))

    }

  }

  dis.mat

}

# Distance based on stationary Ornstein-Uhlenback model
#' @rdname distances
dist.sou = function (meanRPKM = NULL) {

  object_n <- ncol(meanRPKM)
  #gene_n <- nrow(meanRPKM)

  dis.mat <- matrix(0, nr = object_n, nc = object_n)


  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      V11 <- var(log2(meanRPKM[,i]+1))
      V22 <- var(log2(meanRPKM[,j]+1))
      V12 <- cov(log2(meanRPKM[,i]+1), log2(meanRPKM[,j]+1))

      dis.mat[j,i] <- -log(V12/sqrt(V11*V12))

    }

  }

  dis.mat

}

# Converntional expression distance
#' @rdname distances
dist.ced = function (meanRPKM = NULL) {

  object_n <- ncol(meanRPKM)
  gene_n <- nrow(meanRPKM)

  dis.mat <- matrix(0, nr = object_n, nc = object_n)


  for (i in 1:(object_n-1)) {

    for (j in (i+1):object_n) {

      V11 <- var(log2(meanRPKM[,i]+1))
      V22 <- var(log2(meanRPKM[,j]+1))
      V12 <- cov(log2(meanRPKM[,i]+1), log2(meanRPKM[,j]+1))

      dis.mat[j,i] <- V11+V22-2*V12
      #dis.mat[j,i] <- (euc.dist(log2(meanRPKM[,i]+1),log2(meanRPKM[,j]+1)))^2 / gene_n

    }

  }

  dis.mat

}

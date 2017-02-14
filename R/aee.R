#'
#' @title No zero branch length
#'
#' @name no0br
#' @rdname aee
#'
#' @description This function does a small tweak on the tree to remove zero-length branch
#' between root and MRCA of the ingroup
#'
#' @param phy an rooted tree with zero
#'
#' @return returns the tree without zero-length branch
#'
#' @export
no0br = function(phy) {

  if (!inherits(phy, "phylo"))
    stop(paste0(date(),"tree input is not of class \"phylo\""))

  if (is.null(phy$edge.length))
    stop(paste0(date(),": tree has no branch lengths"))

  if (!is.rooted(phy))
    stop(paste0(date(),": tree is not rooted"))

  if (phy$edge.length[1] != 0)
    stop(paste0(date(),": root has no zero-length branch, nothing to do"))

  phy$edge.length[1] <- min(phy$edge.length[-1])
  phy$edge.length[length(phy$edge.length)] <-
    phy$edge.length[length(phy$edge.length)] - phy$edge.length[1]

  phy

}

#'
#'
#' @title Generate an inversed variance matrix from expression profile across species
#'
#' @name varMatInv
#' @rdname aee
#'
#' @description This function generate an inversed variance matrix from expression profiles
#' of one-to-one orthologous genes across species
#'
#' @param objects a vector of objects of class \code{taxonExp} or an object of class \code{taxaExp}
#' @param phy an rooted expression character tree
#' @param taxa oen single character or a vector of characters sepcifying taxa to generate
#' an inversed variance matrix.
#' If one single character "all" is given,
#' all the taxa included in the \code{taxaExp} will be matched and included ("all" by default).
#' @param subtaxa one single character specifying sub taxa to be included in generating
#' an inversed variance matrix
#'
#' @return returns an inversed variance matrix
#'
#' @export
varMatInv = function(objects , phy, taxa = "all", subtaxa) {

  if (!inherits(phy, "phylo"))
    stop(paste0(date(),"tree input is not of class \"phylo\""))

  if (is.null(phy$edge.length))
    stop(paste0(date(),": tree has no branch lengths which is a necessity for \"varMatInv\""))

  if (length(subtaxa) > 1 || subtaxa == "all")
    stop(paste0(date(),": only one subtaxon are allowed here"))

  dismat <- expdist(objects, taxa = taxa, subtaxa = subtaxa, method = "sou") ### using -ln(rho) to estimate pairwise expression distance

  if (!all(row.names(dismat) %in% phy$tip.label ))
    stop(paste0(date(),": taxa or subtaxa names do not match perfectly with tree tip labels, please check them."))

  n_tip <- Ntip(exp_tree)
  n_node <- Nnode(exp_tree)

  if (n_tip != n_node + 1)
    stop(paste0(date(),"tree is not rooted, please make sure tree is properly rooted. "))

  ### extract distances from the tree
  nodes_dist <- dist.nodes(phy)
  corrmat <- apply(nodes_dist, c(1,2), function(x) exp(-x))

  ### stationary variance
  exp_table <- exptabTE(objects, taxa = "all", subtaxa = subtaxa, logrithm = T)
  stat_var <- mean(apply(exp_table, 2, var))
  var_corrmat <- stat_var * corrmat

  solve(var_corrmat)

}

#'
#' @title Ancestral Expression Estimation
#'
#' @name aee
#' @rdname aee
#'
#' @description This function esitmates ancestral expression profile and related statistical uncertainty
#'
#' @param x a vector of known expression profile, preferably log-transformed expression levels (e.g. log RPKM)
#' @param phy a phylogenetic tree in the form of object "phylo"
#' @param mat a matrix generated from "varMatInv" function
#' @param CI a logical specifying whether to return the 95% confidence intervals
#' of the estimated ancestral expression levels
#'
#' @return returns a list containing estimated ancestral expression profile
#' as well as other requested parameters
#'
#' @export
aee = function(x, phy, mat, CI = TRUE) {

  ### checking input formats
  if (!inherits(phy,"phylo"))
    stop(paste0(date(),": \"phy\" input is not of class \"phylo\""))

  if (is.null(phy$edge.length))
    stop(paste0(date(),": tree has no branch lengths which is a necessity for \"aee\""))

  if (!is.null(names(x))) {
    if (all(names(x) %in% phy$tip.label))
      x <- x[phy$tip.label]
    else warining(paste0(date(),
      "characters do not match perfectly between expression profile vector names and tree tip labels,
      only using tree tip labels in the following analysis"))
  }

  ### checking if the tree is rooted
  n_tip <- Ntip(exp_tree)
  n_node <- Nnode(exp_tree)

  if (n_tip != n_node + 1)
    stop(paste0(date(),"tree is not rooted, please make sure tree is properly rooted. "))

  ancestral <- list()

  expr <- numeric(length = n_tip + n_node) ### expression values vector initiate

  expr[1:n_tip] <- if (is.null(names(x))) x else x[phy$tip.label] ### given expression values of tips
  expr[(n_tip+1):(n_tip+n_node)] <- NA ### ancestral nodes expression values to be estimated

  tr_edges <- phy$edge

  ### using children nodes to estimate ancestral nodes' expression

  while (any (is.na(expr[(n_tip+2):(n_tip+n_node)]))) {

    ### looping while ancestral expression values are not computed at all internal nodes except for root

    for (i in (n_tip+2):(n_tip+n_node)) {

      if (is.na(expr[i])) { ### node that has not been computed yet

        child_nodes <- tr_edges[tr_edges[,1] == i, 2]

        if (any(is.na(expr[child_nodes]))) ### bypass empty children nodes
          next

        mu <- mean(expr[child_nodes])

        beta <- unlist(lapply(child_nodes, function(x) - mat[x,i] / mat[i,i]))
        beta0 <- mu * (1 - sum(beta))

        expr[i] <- beta0 + sum(beta * expr[child_nodes])
      }

    }

  }

  ### for the root node

  child_nodes <- tr_edges[tr_edges[,1] == n_tip+1,2]

  mu <- mean(expr[child_nodes])

  beta <- unlist(lapply(child_nodes, function(x) - mat[x,n_tip+1] / mat[n_tip+1,n_tip+1]))
  beta0 <- mu * (1 - sum(beta))

  expr[n_tip+1] <- beta0 + sum(beta * expr[child_nodes])

  ancestral$est <- expr[(n_tip+1):(n_tip+n_node)]

  ### if calculating confidence interval

  if (CI) {

    ci95 <- matrix(nrow = n_node, ncol = 2)

    for (i in (n_tip+1):(n_tip+n_node)) { ### for every node

      tmp = sqrt(1 / mat[i,i]) * qnorm(0.025)
      ci95[(i-n_tip),] = c(expr[i] + tmp, expr[i] - tmp)

    }

    ancestral$ci95 <- ci95
  }

  ancestral
}

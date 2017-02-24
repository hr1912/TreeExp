#'
#' @title Map the expression distance onto a given tree topology
#'
#' @name map.ls
#' @rdname map
#'
#' @param phy a phylogenetic tree
#' @param D a distance matrix has row&column names corresponded to tip labels
#'
#' @description This function computes the least squares brach lengths
#' on a tree with known topology and a distance matrix. Borrowed something from
#' "phytools" package.
#'
#' @return mapped tree and the sum of squares residual error (Q-score)
#'
#' @examples
#' map.tree <- map.ls(phy, D)
#' attr(map.tree, "Q-score")
#'
#' @export
map.ls = function (phy, D) {

    if (!inherits(phy, "phylo"))
        stop(paste0(date(), ": \"phy\" input is not of class \"phylo\""))

    if (!is.binary.tree(phy))
        phy <- multi2di(phy)

    if (is.rooted(phy))
        phy <- unroot(phy)

    ls.tree(phy, D)
}


#'
#' @title No zero branch length
#'
#' @name no0br
#' @rdname no0br
#'
#' @description This function does a small tweak on the tree to remove zero-length branch
#' between root and MRCA of the ingroup, and replace any negative branch lengths
#' with a small length
#'
#' @param phy an rooted tree
#'
#' @return returns the tree without zero-length or negative-length branch
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

    tmp <- sort(phy$edge.length)
    tmp1 <- tmp[tmp>0]

    if (any(tmp < 0)) {

        warning(paste0(date(),": there are negative-length branches in the tree, replacing them with a small length"))
        phy$edge.length[phy$edge.length < 0] <- tmp1[1] / 1e3 ### add a small length for negative-length branches

    }

    phy$edge.length[1] <- tmp1[1]
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

    ### using -ln(rho) to estimate pairwise expression distance
    dismat <- expdist(objects, taxa = taxa, subtaxa = subtaxa, method = "sou")

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
#' @param select indicate if descendents of the node or all tips should be used
#' @param CI a logical specifying whether to return the 95% confidence intervals
#' of the estimated ancestral expression levels
#'
#' @return returns a list containing estimated ancestral expression profile
#' as well as other requested parameters
#'
#' @examples
#'
#' data('tetraexp')
#' dismat <- expdist(tetraexp.objects, taxa = "all", subtaxa = "Brain", method = "sou")
#' exp_tree <- NJ(dismat)
#' exp_tree <- root(exp_tree, outgroup = "Chicken_Brain", resolve.root = T)
#' exp_tree <- no0br(exp_tree)
#' var_mat <- varMatInv(objects = tetraexp.objects,phy = exp_tree,taxa = "all", subtaxa = "Brain")
#' exp_table <- exptabTE(tetraexp.objects, taxa = "all", subtaxa = "Brain")
#' exp_one <- aee(exp_table[1,], exp_tree, var_mat)
#' exp_tree$node.label <- exp_one$est
#' plot(exp_tree)
#'
#' @export
aee = function(x, phy, mat, select = c("all", "descendents") , CI = TRUE) {

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

    expr[1:n_tip] <- if (is.null(names(x))) x else as.numeric(x[phy$tip.label]) ### given expression values of tips
    expr[(n_tip+1):(n_tip+n_node)] <- NA ### ancestral nodes expression values to be estimated

    tr_edges <- phy$edge

    ###

    #select <- match.arg(select)

    ### using children nodes to estimate ancestral nodes' expression

    if (0) { ### direct child nodes problematic? start

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

    } ### problematic? end


    if (select == "descendents") {
    ### estimating ancestral expression with all decendant tips
    for (i in (n_tip+1):(n_tip+n_node)) {

    child_nodes <- getDescendants(phy, i)

    child_tips <- child_nodes[child_nodes < (n_tip+1)] ### only tips

    mu <- mean(expr[child_tips])

    beta <- unlist(lapply(child_tips, function(x) - mat[x,i] / mat[i,i]))
    beta0 <- mu * (1 - sum(beta))

    expr[i] <- beta0 + sum(beta * expr[child_tips])

    }

    } else {

    all_tips <- 1:n_tip

    for (i in (n_tip+1):(n_tip+n_node)) {

        mu <- mean(expr[all_tips])

        beta <- unlist(lapply(all_tips, function(x) - mat[x,i] / mat[i,i]))
        beta0 <- mu * (1 - sum(beta))

        expr[i] <- beta0 + sum(beta * expr[all_tips])

    }

    }

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

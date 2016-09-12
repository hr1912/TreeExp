#' @title Neighbor-joining
#'
#' @name NJ
#'
#' @rdname nj
#'
#' @description The famous neighbor-joining tree estimation function from Saitou and Nei (1987).
#'
#' @examples
#'
#' data(tetraexp)
#' dismat <- expdist(tetraexp.objects, taxa = "all",
#'                  subtaxa = "Brain",
#'                  method = "pea")
#' tr <- root(NJ(dismat), "Chicken_Brain")
#' plot(tr)
#'
#' @export
NJ = function (x) {

  reorder(nj(x), "postorder")

}

#' Consensus class using hierarchical clustering
#'
#' Performs hierarchical clustering on a consensus matrix to obtain consensus
#' class labels.
#'
#' @param x a \code{\link{consensus_matrix}}
#' @param k number of clusters
#' @param method linkage type for hierarchical clustering. Defaults to
#'   "average". See \code{\link[stats]{hclust}} for details.
#' @param names vector of sample labels. Default is \code{NULL}, so that
#'   the names are obtained from the \code{hclust} dendrogram labels.
#' @return cluster assignments for the consensus class
#' @family consensus functions
#' @author Derek Chiu
#' @export
#' @examples
#' set.seed(2)
#' x <- replicate(100, rbinom(100, 4, 0.2))
#' cm <- consensus_matrix(x)
#' consensus_class(cm, k = 3)
consensus_class <- function(x, k, method = "average", names = NULL) {
  tree <- hclust(dist(x), method = method)
  cl <- as.factor(cutree(tree, k))
  if (!is.null(names))
    names(cl) <- names
  else
    names(cl) <- tree$labels
  return(cl)
}

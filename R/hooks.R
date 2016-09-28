#' Hook Functions
#'
#' Custom hook functions for AGglomerative NESting and DIvisive ANAlysis
#' clustering algorithms
#'
#' The hierarchical structure returned using \code{agnes} is equivalent to using
#' \code{\link{hclust}}.
#'
#' @param d distance matrix
#' @param k scalar indicating number of clusters to cut tree into
#' @return clustering assignment from \code{agnes} or \code{diana}
#' @name hooks
#' @author Derek Chiu
#' @import cluster bioDist
#' @export
#' @examples
#' library(cluster)
#' data(votes.repub)
#' d <- dist(votes.repub)
#' agnes_hook(d, k = 2)
#' diana_hook(d, k = 2)
agnes_hook <- function(d, k) {
  tmp <- agnes(d, diss = TRUE)
  a <- cutree(tmp, k)
  return(a)
}

#' @rdname hooks
#' @export
diana_hook <- function(d, k) {
  tmp <- diana(d, diss = TRUE)
  a <- cutree(tmp, k)
  return(a)
}

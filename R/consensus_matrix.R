#' Consensus matrix
#'
#' Returns the (weighted) consensus matrix given a data matrix
#'
#' Given a vector of cluster assignments, we first calculate the connectivity
#' matrix and indicator matrix. A connectivity matrix has a 1 if both samples
#' are in the same cluster, and 0 otherwise. An indicator matrix has a 1 if both
#' samples were selected to be used in a subsample of a consensus clustering
#' algorithm, and 0 otherwise. Summation of connectivity matrices and indicator
#' matrices is performed over different subsamples of the data. The consensus
#' matrix is calculated by dividing the aggregated connectivity matrices by the
#' aggregated indicator matrices.
#'
#' If a meta-consensus matrix is desired, where consensus classes of different
#' clustering algorithms are aggregated, we can construct a weighted
#' meta-consensus matrix using \code{weights}.
#'
#' @param data data matrix has rows as samples, columns as replicates
#' @param weights a vector of weights for each algorithm used in meta-consensus
#'   clustering. Must have \code{length(weights)} equal to \code{ncol(data)}.
#' @return a consensus matrix
#' @note When consensus is calculated over bootstrap samples, not every sample
#'   is used in each replication. Thus, there will be scenarios where two
#'   samples are never chosen together in any bootstrap samples. This typically
#'   happens when the number of replications is small. The coordinate in the
#'   consensus matrix for such pairs of samples is \code{NaN} from a 0 / 0
#'   computation. These entries are coerced to 0.
#' @author Derek Chiu
#' @export
#' @examples
#' set.seed(2)
#' x <- replicate(100, rbinom(100, 4, 0.2))
#' w <- rexp(100)
#' w <- w / sum(w)
#' cm1 <- consensus_matrix(x)
#' cm2 <- consensus_matrix(x, weights = w)
consensus_matrix <- function(data, weights = NULL) {
  data <- as.data.frame(data)
  all.IM <- purrr::map(data, indicator_matrix)
  all.CM <- purrr::map(data, connectivity_matrix)
  sum.IM <- Reduce(`+`, all.IM)
  if (!is.null(weights)) {
    weighted.CM <- purrr::map2(all.CM, weights, `*`)
    sum.CM <- Reduce(`+`, weighted.CM) * length(weights)
  } else {
    sum.CM <- Reduce(`+`, all.CM)
  }
  Reduce(`/`, list(sum.CM, sum.IM)) %>%
    magrittr::inset(is.nan(.), 0)
}

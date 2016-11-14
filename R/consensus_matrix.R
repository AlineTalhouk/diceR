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
#' @author Derek Chiu
#' @export
#' @examples
#' set.seed(2)
#' x <- replicate(100, rbinom(100, 4, 0.2))
#' w <- rexp(100)
#' w <- w / sum(w)
#' consensus_matrix(x)
#' consensus_matrix(x, weights = w)
consensus_matrix <- function(data, weights = NULL) {
  all.IM <- plyr::alply(data, 2, indicator_matrix)
  all.CM <- plyr::alply(data, 2, connectivity_matrix)
  sum.IM <- Reduce('+', all.IM)
  if (!is.null(weights)) {
    weighted.CM <- mapply('*', all.CM, weights, SIMPLIFY = FALSE)
    sum.CM <- Reduce('+', weighted.CM) * length(weights)
  } else {
    sum.CM <- Reduce('+', all.CM)
  }
  cons.mat <- Reduce('/', list(sum.CM, sum.IM))
  return(cons.mat)
}

#' Connectivity matrix
#' @noRd
connectivity_matrix <- function(cls) {
  # cls is a vector of cluster assignments
  cm <- cls %>%
    rep(., length(.)) %>%
    matrix(ncol = sqrt(length(.)))
  for (j in 1:ncol(cm)) {
    if (is.na(cm[j, j])) {
      cm[, j] <- 0
    } else {
      cm[, j] <- ifelse(cm[j, j] != cm[, j] | is.na(cm[, j]), 0, 1)
    }
  }
  rownames(cm) <- colnames(cm) <- names(cls)
  return(cm)
}

#' Indicator matrix
#' @noRd
indicator_matrix <- function(cls) {
  # cls is a vector of cluster assignments
  im <- cls %>%
    rep(., length(.)) %>%
    matrix(ncol = sqrt(length(.)))
  for (j in 1:ncol(im)) {
    if (is.na(im[j, j])) {
      im[, j] <- 0
    } else {
      im[, j] <- ifelse(is.na(im[, j]), 0, 1)
    }
  }
  rownames(im) <- colnames(im) <- names(cls)
  return(im)
}

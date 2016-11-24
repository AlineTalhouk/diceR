#' K-modes
#' 
#' Combine clustering results using K-modes.
#' 
#' Combine clustering results generated using different algorithms and
#' different data perturbations by k-modes. This method is the categorical data
#' analog of k-means clustering. Complete cases are needed: i.e. no \code{NA}s.
#' If the matrix contains \code{NA}s those are imputed by majority voting
#' (after class relabeling).
#' 
#' @param E a matrix of clusterings with number of rows equal to the number of
#' cases to be clustered, number of columns equal to the clustering obtained by
#' different resampling of the data, and the third dimension are the different
#' algorithms. Matrix may already be two-dimensional.
#' @param is.relabelled logical; if \code{FALSE} the data will be relabelled
#' using the first clustering as the reference.
#' @param seed random seed for reproducibility
#' @return a vector of cluster assignments based on k-modes
#' @family consensus functions
#' @author Aline Talhouk
#' @export
#' @examples
#' # Calculate for a fraction of first algorithm
#' data(E_imputed)
#' table(k_modes(E_imputed[1:100, , 1, drop = FALSE], is.relabelled = FALSE))
k_modes <- function(E, is.relabelled = TRUE, seed = 1) {
  set.seed(seed)
  # flatten (and relabel) E
  flat_E <- flatten_E(E, is.relabelled = is.relabelled)
  # fill in missing values if any using majority voting
  if (anyNA(flat_E)) {
    flat_E <- t(apply(flat_E, 1, function(x) {
      x[which(is.na(x))] <- names(which.max(table(x)))
      return(x)
    }))
  }
  # k-modes clustering
  k_modes <- klaR::kmodes(flat_E,
                          modes = max(unlist(apply(flat_E, 2, function(x)
                            length(names(table(x)))))))
  return(k_modes$cluster)
}

#' Majority voting
#' 
#' Combine clustering results using majority voting.
#' 
#' Combine clustering results generated using different algorithms and
#' different data perturbations by majority voting. The class of a sample is
#' the cluster label which was selected most often across algorithms and
#' subsamples.
#' 
#' @param E a matrix of clusterings with number of rows equal to the number of
#' cases to be clustered, number of columns equal to the clustering obtained by
#' different resampling of the data, and the third dimension are the different
#' algorithms. Matrix may already be two-dimensional.
#' @param is.relabelled logical; if \code{FALSE} the data will be relabelled
#' using the first clustering as the reference.
#' @return a vector of cluster assignments based on majority voting
#' @family consensus functions
#' @author Aline Talhouk
#' @export
#' @examples
#' data(E_imputed)
#' table(majority_voting(E_imputed, is.relabelled = FALSE))
majority_voting <- function(E, is.relabelled = TRUE) {
  # flatten (and relabel) E
  flat_E <- flatten_E(E, is.relabelled = is.relabelled)
  
  # majority vote
  maj.vote <- apply(flat_E, 1, function(x)
    as.numeric(names(which.max(table(x)))))
  return(maj.vote)
}

#' Cluster-based Similarity Partitioning Algorithm (CSPA)
#'
#' Performs hierarchical clustering on a stack of consensus matrices to obtain
#' consensus class labels.
#'
#' @param E is an array of clustering results.
#' @param k number of clusters
#' @return cluster assignments for the consensus class
#' @family consensus functions
#' @author Derek Chiu
#' @export
#' @examples
#' data(hgsc)
#' dat <- t(hgsc[, -1])[1:100, 1:50]
#' x <- consensus_cluster(dat, nk = 4, reps = 4, algorithms = c("nmf", "hc",
#' "diana"), nmf.method = "lee", progress = FALSE)
#' CSPA(x, k = 4)
CSPA <- function(E, k) {
  assertthat::assert_that(k %in% dimnames(E)[[4]])
  cl <- consensus_combine(E, element = "matrix") %>% 
    magrittr::extract2(as.character(k)) %>% 
    Reduce("+", .) %>% 
    magrittr::divide_by(dim(E)[3]) %>% 
    stats::dist() %>% 
    hc(k = k)
  return(cl)
}

#' Linkage Clustering Ensemble
#' 
#' Generate a cluster assignment from a CTS, SRS, or ASRS similarity matrix.
#' 
#' @param E is an array of clustering results. An error is thrown if there are 
#'   missing values. \code{\link{impute_missing}} can be used beforehand.
#' @param k requested number of clusters
#' @param dcCTS decay constant for CTS matrix
#' @param dcSRS decay constant for SRS matrix
#' @param dcASRS decay constant for ASRS matrix
#' @param R number of repetitions for SRS matrix
#' @param sim.mat similarity matrix; choices are "cts", "srs", "asrs".
#' @family consensus functions
#' @author Johnson Liu
#' @return a vector containing the cluster assignment from either the CTS, SRS, 
#'   or ASRS similarity matrices
#' @export
#' @examples
#' data(hgsc)
#' dat <- t(hgsc[, -1])[1:100, 1:50]
#' x <- consensus_cluster(dat, nk = 4, reps = 4, algorithms = c("nmf",
#' "hc", "diana"), nmf.method = "lee", progress = FALSE)
#' \dontrun{
#' LCE(E = x, k = 4, sim.mat = "asrs")
#' }
#' 
#' x_imputed <- impute_missing(x, dat, nk = 4)$complete
#' LCE(E = x_imputed, k = 4, sim.mat = "cts")
LCE <- function(E, k, dcCTS = 0.8, dcSRS = 0.8, dcASRS = 0.8, R = 10,
                sim.mat = c("cts", "srs", "asrs")) {
  assertthat::assert_that(is.array(E),
                          dcCTS >= 0 && dcCTS <= 1,
                          dcSRS >= 0 && dcSRS <= 1,
                          dcASRS >= 0 && dcASRS <= 1)
  # Check that the Cluster matrix is complete otherwise return Error
  if (anyNA(E)) stop("'E' must be complete for LCE algorithm.")
  S <- switch(match.arg(sim.mat, c("cts", "asrs", "srs")),
              cts = cts(E = E, dc = dcCTS),
              srs = srs(E = E, dc = dcSRS, R = R),
              asrs = asrs(E = E, dc = dcSRS))
  LCE_cl <- hc(stats::dist(S), k)
  return(LCE_cl)
}
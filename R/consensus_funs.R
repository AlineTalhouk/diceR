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
#' data(hgsc)
#' dat <- hgsc[1:100, 1:50]
#' cc <- consensus_cluster(dat, nk = 4, reps = 6, algorithms = "pam", progress =
#' FALSE)
#' table(k_modes(cc[, , 1, 1, drop = FALSE], is.relabelled = FALSE))
k_modes <- function(E, is.relabelled = TRUE, seed = 1) {
  # Flatten and fill in remaining missing entries with majority voting
  mv <- majority_voting(E, is.relabelled = is.relabelled)
  flat_E <- flatten_E(E, is.relabelled = is.relabelled) %>%
    purrr::array_branch(margin = 1) %>%
    purrr::map2(mv, ~ .x %>%
                  magrittr::inset(is.na(.), .y)) %>%
    do.call(rbind, .) %>%
    as.data.frame()
  # k-modes clustering if there are multiple columns of assignments
  if (ncol(flat_E) > 1) {
    set.seed(seed)
    k_modes <- klaR::kmodes(flat_E,
                            modes = max(purrr::map_int(flat_E,
                                                       dplyr::n_distinct)))
    return(k_modes$cluster)
  } else {
    return(flat_E)
  }
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
#' data(hgsc)
#' dat <- hgsc[1:100, 1:50]
#' cc <- consensus_cluster(dat, nk = 4, reps = 6, algorithms = "pam", progress =
#' FALSE)
#' table(majority_voting(cc[, , 1, 1, drop = FALSE], is.relabelled = FALSE))
majority_voting <- function(E, is.relabelled = TRUE) {
  # Flatten (and relabel) E then find most common element in every row
  E %>%
    flatten_E(is.relabelled = is.relabelled) %>%
    apply(1, function(x) as.numeric(names(which.max(table(x)))))
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
#' dat <- hgsc[1:100, 1:50]
#' x <- consensus_cluster(dat, nk = 4, reps = 4, algorithms = c("hc", "diana"),
#' progress = FALSE)
#' CSPA(x, k = 4)
CSPA <- function(E, k) {
  assertthat::assert_that(k %in% dimnames(E)[[4]])
  consensus_combine(E, element = "matrix") %>%
    magrittr::extract2(as.character(k)) %>%
    Reduce(`+`, .) %>%
    magrittr::divide_by(dim(E)[3]) %>%
    stats::dist() %>%
    hc(k = k)
}

#' Linkage Clustering Ensemble
#'
#' Generate a cluster assignment from a CTS, SRS, or ASRS similarity matrix.
#'
#' @param E is an array of clustering results. An error is thrown if there are
#'   missing values. \code{\link{impute_missing}} can be used beforehand.
#' @param k requested number of clusters
#' @param dc decay constant for CTS, SRS, or ASRS matrix
#' @param R number of repetitions for SRS matrix
#' @param sim.mat similarity matrix; choices are "cts", "srs", "asrs".
#' @family consensus functions
#' @author Johnson Liu
#' @return a vector containing the cluster assignment from either the CTS, SRS,
#'   or ASRS similarity matrices
#' @export
#' @examples
#' data(hgsc)
#' dat <- hgsc[1:100, 1:50]
#' x <- consensus_cluster(dat, nk = 4, reps = 4, algorithms = c("km", "hc",
#' "diana"), progress = FALSE)
#' \dontrun{
#' LCE(E = x, k = 4, sim.mat = "asrs")
#' }
#'
#' x <- apply(x, 2:4, impute_knn, data = dat, seed = 1)
#' x_imputed <- impute_missing(x, dat, nk = 4)
#' LCE(E = x_imputed, k = 4, sim.mat = "cts")
LCE <- function(E, k, dc = 0.8, R = 10, sim.mat = c("cts", "srs", "asrs")) {
  assertthat::assert_that(is.array(E), dc >= 0 && dc <= 1)

  # Check that the Cluster matrix is complete otherwise return Error
  if (anyNA(E)) stop("'E' must be complete for LCE algorithm.")
  S <- switch(match.arg(sim.mat, c("cts", "srs", "asrs")),
              cts = cts(E = E, dc = dc),
              srs = srs(E = E, dc = dc, R = R),
              asrs = asrs(E = E, dc = dc))
  hc(stats::dist(S), k)
}

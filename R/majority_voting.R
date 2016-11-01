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

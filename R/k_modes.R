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
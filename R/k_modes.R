#' Combine clustering results using K modes
#' 
#' Combine clustering results generated using different algorithms and different
#' data perturbations
#' 
#' @param E a matrix of clusterings with number of rows equal to the number of
#'   cases to be clustered, number of columns equal to the clustering obtained
#'   by different resampling of the data, and the third dimension is the
#'   different algorithms. Needs complete cases: No NAs. If the matrix contains
#'   NAs those are imputed by majority voting (after class relabeling). Matrix can be flattened already.
#' @param is.relabelled logical; defaults to TRUE, but if FALSE the data will
#'   be relabelled
#' @param is.flat set to TRUE to indicate that impute matrix is already 2D with rows as cases and columns as clusterings 
#' @param seed random seed for reproducibility
#' @return a vector of cluster assignment based on kmodes
#' @author Aline Talhouk
#' @export
#' @examples
#' # Calculate for a fraction of first algorithm
#' data(E_imputed)
#' table(k_modes(E_imputed[1:100, , 1, drop = FALSE], is.relabelled = FALSE))
k_modes <- function(E, is.relabelled = TRUE,is.flat=TRUE, seed = 1) {
  set.seed(seed)
  # take E imputed and reshape into a flat matrix
  flat_E <- E
  if(is.flat==FALSE){
    dim(flat_E) <- c(dim(E)[1],dim(E)[2]*dim(E)[3])}
 
   if (is.relabelled == FALSE) {
    E_f <- apply(flat_E[, -1], 2, function(x) {
      relabel_class(x, flat_E[, 1])
    })
    flat_E <- cbind(flat_E[, 1], E_f)
  }
 
   # Fill in missing values if any using majority voting
  if (anyNA(flat_E)) {
    flat_E <- t(apply(flat_E, 1, function(x) {
      x[which(is.na(x))] <- names(which.max(table(x)))
      return(x)
    }))
  }
  k_modes <- klaR::kmodes(
    flat_E,
    modes = max(unlist(apply(flat_E, 2, function(x)
      length(names(table(x))))))
  )
  return(k_modes$cluster)
}

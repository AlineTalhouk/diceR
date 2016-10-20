#' Combine clustering results using majority voting
#'
#' Combine clustering results generated using different algorithms and different data perturbations
#' @param E a matrix of clusterings with number of rows equal to the number of cases to be clustered, number of columns equal to the clustering obtained by different resampling of the data, and the third dimension is the different algorithms. Matrix can be flattened already
#' @param is.relabelled TRUE/FALSE defaults to TRUE, but if FALSE the data will be relabelled
#' @param is.flat set to TRUE to indicate that impute matrix is already 2D with rows as cases and columns as clusterings 
#' @return a vector of cluster assignment based on majority vote
#' @author Aline Talhouk
#' @export
#' @examples
#' data(E_imputed)
#' table(majority_voting(E_imputed, is.relabelled=FALSE))

majority_voting <- function(E, is.relabelled=TRUE, is.flat=TRUE){
#take E imputed and reshape into a flat matrix
  flat_E <- E
if(is.flat==FALSE){
dim(flat_E) <- c(dim(E)[1],dim(E)[2]*dim(E)[3])}

if(is.relabelled==FALSE){
  E_f <- apply(flat_E[,-1],2,function(x){relabel_class(x,flat_E[,1])})
  flat_E <- cbind(flat_E[,1],E_f)
  }
#majority vote
maj.vote <- as.vector(apply(flat_E, 1, function(x) names(which.max(table(x)))))

return(maj.vote)
}

#' Impute missing values due to subsampling
#'
#' Uses k nearest neighbour to impute missing values due to sampling, not all cases will be imputed using this method. If imputeALL is set to true (default), then cases will be relabelled and missing cases will be called by majority vote
#'
#' @param E a matrix of clusterings with number of rows equal to the number of cases to be clustered, number of columns equal to the clustering obtained by different resampling of the data, and the third dimension is the different algorithms.
#' @param data is data matrix with samples as rows and genes/features as columns
#' @param imputeALL if FALSE the function will only call knn_impute to impute some NAs, if TRUE, the function will use knn_impute and majority voting to eliminate NAs and relabel class
#' @return a flat matrix of clusterings, fully imputed and relabelled.
#' @author Aline Talhouk
#' @export
#' @examples
#' data(hgsc)
#' data <- t(hgsc[,-1])
#' E <- ConClust(data, k = 4, reps = 10, method = c("hcAEucl","kmEucl","scRbf"))
#' imputeMissing(E,data)

imputeMissing <- function(E, data, imputeALL = TRUE) {
  assertthat::assert_that(is.array(E))
  assertthat::assert_that(is.matrix(data))
  assertthat::assert_that(dim(E)[1]==nrow(data))
  assertthat::assert_that(imputeALL==TRUE || imputeALL==FALSE)
  # Flatten the matrix
  E_flat <- E
  dim(E_flat) <- c(dim(E)[1], dim(E)[2] * dim(E)[3])
  # knn impute
  E_imputed <- apply(E_flat, 2, knn_impute, data = data)
  # Relabel and Majority vote
  if (imputeALL == TRUE) {
    E_relabeled <- cbind(E_imputed[, 1],
                         apply(E_imputed[, -1], 2, function(x) {
                           relabel_class(x, E_imputed[, 1])
                         }))
    E_imputed2 <- t(apply(E_relabeled, 1, function(x) {
      x[which(is.na(x))] <- names(which.max(table(x)))
      return(x)
    }))
    return(apply(E_imputed2, 2, as.numeric))
  } else{
    return(E_imputed)
  }
}
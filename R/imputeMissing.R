#' Impute missing values due to subsampling
#' 
#' Uses k nearest neighbour to impute missing values due to sampling. Not all
#' cases will be imputed using this method.
#' 
#' @param E4 is an array of clusterings with number of rows equal to the number of cases to be clustered, number of columns equal to the clustering obtained by different resampling of the data, and the third dimension is the different algorithms and the fourth dimension is the cluster size.
#' @param data is data matrix with samples as rows and genes/features as columns
#' @param imputeALL if \code{FALSE} the function will only call knn_impute to
#'   impute some NAs; if \code{TRUE}, the function will use knn_impute and
#'   majority voting to eliminate NAs and relabel class
#' @return a flat matrix of clusterings, fully imputed and relabelled.
#' @author Aline Talhouk
#' @export
#' @examples
#' data(hgsc)
#' data <- t(hgsc[,-1])[1:200, 1:100]
#' E <- ConClust(data, nc = 2:4, reps = 10, method = c("hcAEucl", "kmEucl",
#' "scRbf"))
#' imputeMissing(E, data)
imputeMissing <- function(E4, data, imputeALL = TRUE) {
  assertthat::assert_that(is.array(E4), is.matrix(data),
                          dim(E4)[1] == nrow(data),
                          imputeALL == TRUE || imputeALL == FALSE)
  # For cluster size, do the following
  E_return <- E4
  for (k in (1:dim(E4)[4])){
  # Flatten the matrix
  E <- E4[,,,k]
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
    E_return[,,,k] <- apply(E_imputed2, 2, as.numeric)
  } else{
    E_return[,,,k] <- E_imputed
  }
  }
  return(E_return)
}
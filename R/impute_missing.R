#' Impute missing values
#' 
#' Impute missing values from bootstrapped subsampling
#' 
#' The default output from \code{consensus_cluster} will undoubtedly contain \code{NA} 
#' entries because each replicate chooses a random subset (with replacement) of 
#' all samples. Missing values are first imputed using KNN (K-Nearest 
#' Neighbours), with the non-missing cases indicating the training set, and
#' missing cases indicating the test set. Not all missing values are guaranteed
#' to be imputed by KNN. See \code{\link[class]{knn}} for details. Thus, any
#' remaining missing values are imputed using majority voting.
#' 
#' @param E 4D array of clusterings from \code{consensus_cluster}. The number of rows is 
#'   equal to the number of cases to be clustered, number of columns is equal to
#'   the clusterings obtained by different resamplings of the data, the third 
#'   dimension are the different algorithms and the fourth dimension are cluster
#'   sizes.
#' @param data data matrix with samples as rows and genes/features as columns
#' @param nk cluster size to extract data for (single value)
#' @param seed random seed for KNN reproducibility
#' @return A list with the following elements
#'   \item{knn}{E with (some) missing cases imputed using KNN}
#'   \item{complete}{flattened matrix of clusterings with complete cases imputed
#'   using KNN and majority voting, and relabelled, for chosen \code{k}}
#' @note We consider 5 nearest neighbours and the minimum vote for definite 
#'   decision is 3.
#' @author Aline Talhouk
#' @export
#' @examples
#' data(hgsc)
#' data <- t(hgsc[, -1])[1:100, 1:50]
#' E <- consensus_cluster(data, nk = 3:4, reps = 10, algorithms = c("hc", "km",
#' "sc"), progress = FALSE)
#' sum(is.na(E))
#' E_imputed <- impute_missing(E, data, 4)
#' sum(is.na(E_imputed))
impute_missing <- function(E, data, nk, seed = 123456) {
  assertthat::assert_that(is.array(E), is.matrix(data),
                          dim(E)[1] == nrow(data),
                          as.character(nk) %in% dimnames(E)[[4]])
  idk <- match(nk, dimnames(E)[[4]])
  # Flatten the matrix
  E_relabeled <- flatten_E(E[, , , idk, drop = FALSE], is.relabelled = FALSE) %>% 
    extract(, apply(., 2, function(x) !all(is.na(x))), drop = FALSE)
  # Relabel and Majority vote
  if (ncol(E_relabeled) > 1) {
    E_complete <- t(apply(E_relabeled, 1, function(x) {
      x[which(is.na(x))] <- as.numeric(names(which.max(table(x))))
      return(x)
    }))
    return(E_complete)
  } else {
    return(E_relabeled)  
  }
}

#' K-Nearest Neighbours imputation
#' Set NA cases as test, otherwise training; neighbours = 5, min vote = 3
#' @param x clustering object
#' @param data data matrix
#' @param seed random seed for knn imputation
#' @export
impute_knn <- function(x, data, seed) {
  set.seed(seed)
  ind <- !is.na(x)
  if (any(ind))
    x[!ind] <- class::knn(train = data[ind, ], test = data[!ind, ], cl = x[ind],
                          k = 5, l = 3)
  return(x)
}
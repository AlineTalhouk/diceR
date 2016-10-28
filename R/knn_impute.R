#' K-Nearest Neighbours Imputation
#'
#' Use K-Nearest Neighbours to fill in samples not clustered due to bootstrap
#' resampling in consensus clustering.
#'
#' The number of neighbours considered is 5 (default) and the minimum vote for
#' definite decision is 3 (default).
#'
#' @param res a vector representing a column of the output from \code{ConClust}
#' @param data data matrix with samples as rows and genes as columns
#' @param seed random seed for reproducibility
#' @return results matrix with missing values filled in
#' @author Aline Talhouk
#' @export
#' @examples
#' # Generate data from ConClust
#' set.seed(2)
#' x <- replicate(100, rnorm(100))
#' CC <- ConClust(x, nc = 4, reps = 10, method = "scRbf")
#' orig.cl <- CC[, 1, 1, 1]
#' sum(is.na(orig.cl))
#'
#' # Still missing values because there is doubt in the voting for some cases
#' impute.cl <- knn_impute(orig.cl, x)
#' sum(is.na(impute.cl))
knn_impute <- function(res, data, seed = 123456) {
  set.seed(seed)
  res.n <- res
  ind <- is.na(res)
  cl <- res[!ind]
  train <- data[!ind, ]
  test <- data[ind, ]
  knn.res <- class::knn(train, test, cl, k = 5, l = 3, prob = TRUE)
  res.n[ind] <- knn.res
  return(res.n)
}

#' Impute missing values
#'
#' Impute missing values from bootstrapped subsampling
#'
#' The default output from `consensus_cluster` will undoubtedly contain `NA`
#' entries because each replicate chooses a random subset (with replacement) of
#' all samples. Missing values should first be imputed using [impute_knn()]. Not
#' all missing values are guaranteed to be imputed by KNN. See [class::knn()]
#' for details. Thus, any remaining missing values are imputed using majority
#' voting.
#'
#' @param E 4D array of clusterings from `consensus_cluster`. The number of rows
#'   is equal to the number of cases to be clustered, number of columns is equal
#'   to the clusterings obtained by different resamplings of the data, the third
#'   dimension are the different algorithms and the fourth dimension are cluster
#'   sizes.
#' @param data data matrix with samples as rows and genes/features as columns
#' @param nk cluster size to extract data for (single value)
#' @return If flattened matrix consists of more than one repetition, i.e. it
#'   isn't a column vector, then the function returns a matrix of clusterings
#'   with complete cases imputed using majority voting, and relabelled, for
#'   chosen `k`.
#' @author Aline Talhouk
#' @family imputation functions
#' @export
#' @examples
#' data(hgsc)
#' dat <- hgsc[1:100, 1:50]
#' E <- consensus_cluster(dat, nk = 3:4, reps = 10, algorithms = c("hc", "km",
#' "sc"), progress = FALSE)
#' sum(is.na(E))
#' E_imputed <- impute_missing(E, dat, 4)
#' sum(is.na(E_imputed))
impute_missing <- function(E, data, nk) {
  assertthat::assert_that(is.array(E),
                          dim(E)[1] == nrow(data),
                          as.character(nk) %in% dimnames(E)[[4]])
  idk <- match(nk, dimnames(E)[[4]])
  # Flatten the matrix
  E_relabeled <- flatten_E(E[, , , idk, drop = FALSE],
                           is.relabelled = FALSE) %>%
    magrittr::extract(, apply(., 2, function(x) !all(is.na(x))), drop = FALSE)
  # Relabel and Majority vote
  if (ncol(E_relabeled) > 1) {
    E_complete <- apply(E_relabeled, 1, function(x) {
      x[which(is.na(x))] <- as.numeric(names(which.max(table(x))))
      return(x)
    }) %>% t()
    return(E_complete)
  } else {
    return(E_relabeled)
  }
}

#' K-Nearest Neighbours imputation
#'
#' The non-missing cases indicate the training set, and missing cases indicate
#' the test set.
#'
#' @param x clustering object
#' @param data data matrix
#' @param seed random seed for knn imputation reproducibility
#' @return An object with (potentially not all) missing values imputed with
#'   K-Nearest Neighbours.
#' @author Aline Talhouk
#' @family imputation functions
#' @note We consider 5 nearest neighbours and the minimum vote for definite
#'   decision is 3.
#' @export
#' @examples
#' data(hgsc)
#' dat <- hgsc[1:100, 1:50]
#' x <- consensus_cluster(dat, nk = 4, reps = 4, algorithms = c("km", "hc",
#' "diana"), progress = FALSE)
#' x <- apply(x, 2:4, impute_knn, data = dat, seed = 1)
impute_knn <- function(x, data, seed = 123456) {
  set.seed(seed)
  ind <- !is.na(x)
  if (any(ind))
    x[!ind] <- class::knn(train = data[ind, ], test = data[!ind, ], cl = x[ind],
                          k = 5, l = 3)
  return(x)
}

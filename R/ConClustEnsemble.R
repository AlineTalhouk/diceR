#' Function to perform clustering using \code{\link{ConClust}} with all default 
#' methods using user-sepcified number of clusters and number of repetitions.
#' 
#' @param X data set with rows as samples, columns as genes
#' @param k number of clusters in the base clusterings
#' @param reps number of repetitions
#' @param method vector of clustering algorithms for performing consensus 
#'   clustering. Must be any number of the following: "nmfDiv", "nmfEucl", 
#'   "hcAEucl", "hcDianaEucl", "kmEucl", "kmSpear", "pamEucl", "pamSpear", 
#'   "apEucl", "scRbf", "gmmBIC", "biclust".
#'   
#' @return matrix of cluster ensemble
#' @author Johnson Liu
#' @importFrom stats complete.cases
#' @export
#' @examples 
#' # Generate a 489 by 10 clustering ensemble using 2 algorithms in ConClust
#' data(hgsc)
#' dat <- hgsc[, -1]
#' rownames(dat) <- hgsc[, 1]
#' E <- ConClustEnsemble(X = t(dat), k = 4, reps = 5, method = c("hcAEucl", "kmEucl"))
ConClustEnsemble <- function(X, k, reps = 1000, method = NULL) {
  assertthat::assert_that(is.matrix(X))
  assertthat::assert_that(reps > 1)  # can't compute consensus matrix on one rep
  df_conClust <- data.frame(t(X))
  colnames(df_conClust) <- rownames(X)
  rownames(df_conClust) <- colnames(X)
  method.list <- c("nmfDiv", "nmfEucl", "hcAEucl", "hcDianaEucl", "kmEucl",
                   "kmSpear", "pamEucl", "pamSpear", "apEucl", "scRbf",
                   "gmmBIC", "biclust")
  if (is.null(method)) {
    method <- method.list
  } else {
    method.input <- method
    method <- match.arg(method.input, method.list, several.ok = TRUE)
    if (length(method) < length(method.input))
      stop("Not all elements in 'method' have a partial match.")
  }
  E <- ConClust(df_conClust, k = k, reps = reps, method = method)
  if (is.na(dim(E)[3])) {
    E_imputed <- array(dim = c(dim(E), 1))
  } else{
    E_imputed <- array(dim = dim(E))
  }
  for (i in 1:dim(E_imputed)[3]) {
    E_imputed[, , i] <- apply(E[, , i], 2, knn_impute, data = X)
  }
  E_cs <- consensus_summary(E_imputed, k = k)
  for (i in 2:length(E_cs)) {
    E_cs[[i]]$consensus_class <- relabel_class(E_cs[[1]]$consensus_class,
                                               E_cs[[i]]$consensus_class)
  }
  for (i in 1:length(E_cs)) {
    E_imputed[, , i] <-
      apply(E_imputed[, , i], 2, relabel_class,
            cl.ref = E_cs[[i]]$consensus_class)
  }
  E_imputed_relabelled <- E_imputed
  flat_E <- E_imputed_relabelled
  # E_imputed_relabelled_voted <- array(dim = dim(E_imputed_relabelled))
  # #If there are still NAs after imputing with KNN, then do majority vote
  # if (!(sum(!complete.cases(E_imputed_relabelled)) == 0)) {
  #   E_imputed_relabelled_voted <- majVote(E_imputed_relabelled)
  # } else{
  #   E_imputed_relabelled_voted <-
  #     apply(E_imputed_relabelled_voted, c(1, 2, 3), as.numeric)
  # }
  dim(flat_E) <- c(dim(E_imputed_relabelled)[1],
                   dim(E_imputed_relabelled)[2] * dim(E_imputed_relabelled)[3])
  if (anyNA(flat_E)) {
    flat_E <- t(apply(flat_E, 1, function(x) {
      x[which(is.na(x))] <- names(which.max(table(x)))
      return(x)
    }))
  }
  flat_E<-apply(flat_E,c(1,2),as.numeric)
  return(flat_E)
}

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

  E_imputed <- array(dim = dim(E))
  
  #Use k nearest neighbour to get rid of some NAs
  for (i in 1:dim(E_imputed)[3]) {
    E_imputed[, , i] <- apply(E[, , i], 2, knn_impute, data = X)
  }
  #Get consensus summary
  E_cs <- consensus_summary(E_imputed, k = k)
  if(length(E_cs)>1){
    for (i in 2:length(E_cs)) {
      E_cs[[i]]$consensus_class <- relabel_class(E_cs[[1]]$consensus_class,
                                                 E_cs[[i]]$consensus_class)
    }
  }
  #Relabel and flatten the layers
  for (i in 1:length(E_cs)) {
    E_imputed[, , i] <-
      apply(E_imputed[, , i], 2, relabel_class,
            cl.ref = E_cs[[i]]$consensus_class)
  }
  E_imputed_relabelled <- E_imputed
  flat_E <- E_imputed_relabelled
  dim(flat_E) <- c(dim(E_imputed_relabelled)[1],
                   dim(E_imputed_relabelled)[2] * dim(E_imputed_relabelled)[3])
  #Get rid of remaining NAs with majority voting
  if (anyNA(flat_E)) {
    flat_E <- t(apply(flat_E, 1, function(x) {
      x[which(is.na(x))] <- names(which.max(table(x)))
      return(x)
    }))
  }
  flat_E<-apply(flat_E,c(1,2),as.numeric)
  return(flat_E)
}

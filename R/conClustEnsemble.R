#' Title : function to generate a cluster ensemble from the data set X using K-means algorithm
#'
#' @param X : a data set, rows are observations, columns are variables
#' @param M : number of base clusterings in the ensemble
#' @param k : number of clusters in the base clusterings
#' @param scheme : cluster ensemble generating scheme (1=fixed k, 2=random k)
#'
#' @return E: matrix of cluster ensemble
#' @export
conClustEnsemble <- function(X, k, reps = 1000) {
  assertthat::assert_that(is.matrix(X))
  E <- NULL
  df_conClust <- data.frame(t(X))
  colnames(df_conClust) <- rownames(X)
  rownames(df_conClust) <- colnames(X)
  E <-
    ConClust::ConClust(
      df_conClust,
      k = k,
      reps = reps,
      method = c(
        "nmfDiv",
        "nmfEucl",
        "hcAEucl",
        "hcDianaEucl",
        "kmEucl",
        "kmSpear",
        "pamEucl",
        "pamSpear",
        "apEucl",
        "scRbf",
        "gmmBIC",
        "biclust"
      )
    )
  if (is.na(dim(E)[3])) {
    E_imputed <- array(dim = c(dim(E), 1))
  } else{
    E_imputed <- array(dim = dim(E))
  }
  for (i in 1:dim(E_imputed)[3]) {
    E_imputed[, , i] <- apply(E[, , i], 2, knnImpute, tdat = X)
  }
  E_cs <- consensusSummaries(E_imputed, k = k)
  for(i in 2:length(E_cs)){
    E_cs[[i]]$consensusClass<-classRelabel(E_cs[[1]]$consensusClass,E_cs[[i]]$consensusClass)
  }
  for (i in 1:length(E_cs)) {
    E_imputed[, , i] <-
      apply(E_imputed[, , i], 2, classRelabel, cl.ref = E_cs[[i]]$consensusClass)
  }
  E_imputed_relabelled <- E_imputed
  E_imputed_relabelled_voted <- array(dim = dim(E_imputed_relabelled))
  #If there are still NAs after imputing with KNN, then do majority vote
  if (!(sum(!complete.cases(E_imputed_relabelled)) == 0)) {
    E_imputed_relabelled_voted <- majVote(E_imputed_relabelled)
  } else{
    E_imputed_relabelled_voted <-
      apply(E_imputed_relabelled_voted, c(1, 2, 3), as.numeric)
  }
  return(E_imputed_relabelled_voted)
}

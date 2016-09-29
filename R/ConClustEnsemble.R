#' Generate a cluster ensemble using the K-means algorithm
#'
#' Function to perform clustering using ConClust::ConClust with default methods nmfDiv,nmfEucl,hcAEucl,hcDianaEucl, kmEucl,kmSpear,pamEucl,apEucl,
#' apEucl, scRbf, gmmBIC, and biclust with user-sepcified number of clusters and number of repetitions
#'
#' @param X a data set, rows are observations, columns are variables
#' @param M number of base clusterings in the ensemble
#' @param k number of clusters in the base clusterings
#'
#' @return matrix of cluster ensemble
#' @author Johnson Liu
#' @importFrom stats complete.cases
#' @export
#' @examples 
#' library(magrittr)
#' library(tidyr)
#' library(reshape2)
#' library(biostatUtil)
#' library(ConClust)
#' data("hgsc")
#' dat <- hgsc[,-1]
#' rownames(dat) <- hgsc[,1]
#' tdat <- t(dat)
#' mdat <- melt(tdat) %>%
#'  set_colnames(c("patientID","Gene","Expr"))
#' E<-ConClustEnsemble(X=tdat,k=4,reps=10) 
#' The example will generate a 489 by 120 clustering ensemble using the default 12 algorithms in ConClust
ConClustEnsemble <- function(X, k, reps = 1000) {
  assertthat::assert_that(is.matrix(X))
  E <- NULL
  df_conClust <- data.frame(t(X))
  colnames(df_conClust) <- rownames(X)
  rownames(df_conClust) <- colnames(X)
  E <-
    ConClust(
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
    E_imputed[, , i] <- apply(E[, , i], 2, knn_impute, data = X)
  }
  E_cs <- consensus_summary(E_imputed, k = k)
  for(i in 2:length(E_cs)){
    E_cs[[i]]$consensus_class<-relabel_class(E_cs[[1]]$consensus_class,E_cs[[i]]$consensus_class)
  }
  for (i in 1:length(E_cs)) {
    E_imputed[, , i] <-
      apply(E_imputed[, , i], 2, relabel_class, cl.ref = E_cs[[i]]$consensus_class)
  }
  E_imputed_relabelled<-E_imputed
  flat_E <- E_imputed_relabelled
  # E_imputed_relabelled_voted <- array(dim = dim(E_imputed_relabelled))
  # #If there are still NAs after imputing with KNN, then do majority vote
  # if (!(sum(!complete.cases(E_imputed_relabelled)) == 0)) {
  #   E_imputed_relabelled_voted <- majVote(E_imputed_relabelled)
  # } else{
  #   E_imputed_relabelled_voted <-
  #     apply(E_imputed_relabelled_voted, c(1, 2, 3), as.numeric)
  # }
  dim(flat_E) <- c(dim(E_imputed_relabelled)[1],dim(E_imputed_relabelled)[2]*dim(E_imputed_relabelled)[3])
  if(anyNA(flat_E)){
    flat_E <- t(apply(flat_E,1,
                      function(x){x[which(is.na(x))] <- names(which.max(table(x)));x}))
  }
  class(flat_E)<-"numeric"
  return(flat_E)
}

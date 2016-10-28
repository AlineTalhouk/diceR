#' Runs through the full diceR algorithm
#'
#' A function that executes the diceR algorithm
#'
#' (Needs to be more informative) A function that runs through the full diceR algorithm
#'
#' @param data a data set with rows as observations, columns as variables
#' @param nk a rangeof cluster sizes (or a single value).
#' @param R number of data subsamples to obtain. (See ConClust documentation)
#' @param algorithms clustering algorithms to be used in the ensemble. Current options are c("nmfDiv", "nmfEucl", "hcAEucl", "hcDianaEucl", "kmEucl","kmSpear", "pamEucl", "pamSpear", "apEucl", "scRbf", "gmmBIC","biclust"). Please see ConClust for definition.
#' @param consensusFUNS Consensus Functions to use. Current options are majority voting, k-modes, CSPA, LCE
#' @param prune a logical value to indicate whether algorithm pruning is needed prior to consensus clustering. defaults to FALSE.
#' @param weigh a logical value to indicate whether after pruning, certain algorithms should have more weight than others. defaults to FALSE.If weigh is TRUE pruning will be set to TRUE as well.
#' @param evaluate indicate whether internal evaluation is needed. evaluate=c("internal","external","both"). If evaluate is either external or both, then a reference class must be provided.
#' @param refClass a vector of length n, indicating the reference class to compare against.
#' @export
dice <- function(data,
                 nk,
                 R = 10,
                 algorithms = c("hcAEucl", "kmEucl", "scRbf", "gmmBIC"),
                 consensusFUNS = c("kmodes", "CSPA", "majority", "LCE"),
                 prune = FALSE,
                 weigh = FALSE,
                 evaluate = "internal",
                 refClass = NULL) {
  # Check that inputs are correct
  if (length(dim(data)) != 2) {
    stop("Data should be two dimensional")
  }
  if (prune == FALSE & weigh == TRUE) {
    prune <- TRUE
  }
  if ((evaluate == "external" |
       evaluate == "both") &
      !is.null(refClass)) {
    stop("Reference Class should be imputted or set evaluate to 'internal'")
  }
  
  n <- dim(data)[1]
  ncf <- length(consensusFUNS)
  
  # Generate Diverse Cluster Ensemble
  E <- ConClust(data,
                nc = nk,
                reps = R,
                method = algorithms)
  
  # Pruning
  if (prune == TRUE) {
    #Function that Derek writing
    #Must return Enew with nalgs< for original E
  if (reweigh == TRUE) {
    #Future function
    #Must return Enew with nalgs< for original E and certain slices are assigned more weight
  }
  } else {
    Enew <- E
  }
  
  # Impute Missing Values using KNN and majority vote
  Ecomp <- imputeMissing(Enew, data, imputeALL = TRUE)

  if (!is.null(refClass)) {
    
  }
  Final <- matrix(NA, nrow = n, ncol = ncf)
  for (i in 1:ncf) {
    Final[, i] <- switch(
      consensusFUNS[i],
      kmodes = k_modes(Ecomp$E_imputed2),
      majority = majority_voting(Ecomp$E_imputed2),
      CSPA = majority_voting(Ecomp$E_imputed2),  # place holder
      LCE = LCE(drop(Ecomp$E_imputed2), nk)
    )
  }
  
  # Relabel Final Clustering
  if (ncf == 1) {
    # if one consensus algorithm is requested only no need to relabel
    FinalR <- Final
  } else if (ncf == 2) {
    FinalR <-
      cbind(Final[, 1], as.numeric(relabel_class(Final[, 2], Final[, 1])))
  } else{
    FinalR <- cbind(Final[, 1],
        apply(Final[,-1], 2, function(x) {
          as.numeric(relabel_class(x, Final[, 1]))
        }))
  }
  colnames(FinalR) <- consensusFUNS
  rownames(FinalR) <- rownames(data)
  return(list(clusters = FinalR))
}
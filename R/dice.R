#' Diverse Clustering Ensemble
#' 
#' Runs consensus clustering across subsamples, algorithms, and number of 
#' clusters (k).
#' 
#' @param data matrix with rows as observations, columns as variables
#' @param nk a range of cluster sizes (or a single value)
#' @param R number of data subsamples to generate. See \code{\link{ConClust}} 
#'   for details.
#' @param algorithms clustering algorithms to be used in the ensemble. Current 
#'   options are "nmfDiv", "nmfEucl", "hcAEucl", "hcDianaEucl", "kmEucl", 
#'   "kmSpear", "pamEucl", "pamSpear", "apEucl", "scRbf", "gmmBIC", and 
#'   "biclust". See \code{\link{ConClust}} for details.
#' @param consensusFUNS consensus functions to use. Current options are "kmodes"
#'   (k-modes), "majority" (majority voting), "CSPA" (hierarchical clustering), 
#'   "LCE" (linkage clustering ensemble)
#' @param sim.mat type of similarity matrix. One of "cts", "srs", "asrs. See
#'   \code{\link{LCE}} for details.
#' @param trim logical; if \code{TRUE}, the number of algorithms in 
#'   \code{algorithms} is reduced based on internal validity index performance 
#'   prior to consensus clustering by \code{consensusFUNS}. Defaults to 
#'   \code{FALSE}.
#' @param reweigh logical; if \code{TRUE}, algorithms are reweighted based on 
#'   internal validity index performance after trimming. Well-performing 
#'   algorithms are given higher weight prior to consensus clustering by 
#'   \code{consensusFUNS}. Defaults to \code{FALSE}. Ignored if \code{trim = 
#'   FALSE}.
#' @param evaluate logical; if \code{TRUE} (default), validity indices are 
#'   returned. Internal validity indices are always computed. If \code{refClass}
#'   is not \code{NULL}, then external validity indices will also be computed.
#' @param refClass reference class; a vector of length equal to the number of 
#'   observations.
#' @param progress logical; if \code{TRUE} (default), progress bar is shown.
#' @return A final clustering assignment from the diverse clustering ensemble
#'   method.
#' @author Aline Talhouk, Derek Chiu
#' @export
#' @examples 
#' library(dplyr)
#' data(hgsc)
#' dat <- t(hgsc[, -1])
#' refClass <- data.frame(initCol = rownames(dat)) %>%
#' tidyr::separate(initCol,
#'                 into = c("patientID", "Class"),
#'                 sep = "_") %>% 
#'   magrittr::use_series(Class) %>% 
#'   factor() %>% 
#'   as.integer()
#' dice.obj <- dice(dat, nk = 4, algorithms = c("hcAEucl", "hcDianaEucl",
#' "pamEucl", "pamSpear"), consensusFUNS = c("kmodes", "majority", "LCE"), trim
#' = TRUE, refClass = refClass)
#' str(dice.obj)
dice <- function(data, nk, R = 10,
                 algorithms = c("hcAEucl", "kmEucl", "scRbf", "gmmBIC"),
                 consensusFUNS = c("kmodes", "CSPA", "majority", "LCE"),
                 sim.mat = c("cts", "srs", "asrs"),
                 trim = FALSE, reweigh = FALSE, evaluate = TRUE,
                 refClass = NULL, progress = TRUE) {
  
  # Check that inputs are correct
  assertthat::assert_that(length(dim(data)) == 2)
  if (!is.null(refClass)) assertthat::assert_that(is.integer(refClass))
  n <- dim(data)[1]
  ncf <- length(consensusFUNS)
  
  # Generate Diverse Cluster Ensemble
  E <- ConClust(data, nc = nk, reps = R, method = algorithms,
                progress = progress)
  
  # Evaluate, trim, and reweigh
  res.obj <- consensus_trim(data, E, ref.cl = refClass, reweigh = reweigh)
  if (trim) E <- res.obj$data.new
  
  # Select k
  k <- res.obj$eval$k
  
  # Impute Missing Values using KNN and majority vote
  imp.obj <- impute_missing(E, data)
  Ecomp <- imp.obj$complete
  
  # Consensus functions
  Final <- matrix(NA, nrow = n, ncol = ncf,
                  dimnames = list(rownames(data), consensusFUNS))
  for (i in 1:ncf) {
    Final[, i] <- switch(consensusFUNS[i],
                         kmodes = k_modes(Ecomp),
                         majority = majority_voting(Ecomp),
                         CSPA = consensus_class(Ecomp, k), 
                         LCE = LCE(drop(Ecomp), k = k,
                                   sim.mat = match.arg(sim.mat))
    )
  }

  # Relabel Final Clustering using reference
  if (ncf == 1 & is.null(refClass)) {
    # Don't relabel if only one consensus function and no reference class
    FinalR <- Final
  } else {
    # Final classes need to be integer for certain functions to work
    FinalR <- apply(Final, 2, function(x) as.integer(
      relabel_class(cl.pred = x, cl.ref = refClass)))
  }
  
  # Return evaluation output including consensus function results
  if (evaluate)
    eval.obj <- consensus_evaluate(data, E, cons.cl = FinalR,
                                   ref.cl = refClass, plot = FALSE)
  
  # Add the reference class as the first column if provided
  if (!is.null(refClass)) {
    FinalR <- cbind(refClass, FinalR)
  }
  
  return(list(clusters = FinalR, indices = eval.obj))
}
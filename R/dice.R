#' Diverse Clustering Ensemble
#' 
#' Runs consensus clustering across subsamples, algorithms, and number of 
#' clusters (k).
#' 
#' @param data matrix with rows as observations, columns as variables
#' @param nk number of clusters (k) requested; can specify a single integer or a
#'   range of integers to compute multiple k
#' @param reps number of data subsamples to generate. See 
#'   \code{\link{consensus_cluster}} for details.
#' @param algorithms clustering algorithms to be used in the ensemble. Current 
#'   options are "nmfDiv", "nmfEucl", "hcAEucl", "hcDianaEucl", "kmEucl", 
#'   "kmSpear", "pamEucl", "pamSpear", "apEucl", "scRbf", "gmmBIC", and 
#'   "biclust". See \code{\link{consensus_cluster}} for details.
#' @param cons.funs consensus functions to use. Current options are "kmodes" 
#'   (k-modes), "majority" (majority voting), "CSPA" (Cluster-based Similarity 
#'   Partitioning Algorithm), "LCE" (linkage clustering ensemble)
#' @param sim.mat type of similarity matrix. One of "cts", "srs", "asrs. See 
#'   \code{\link{LCE}} for details.
#' @param trim logical; if \code{TRUE}, the number of algorithms in 
#'   \code{algorithms} is reduced based on internal validity index performance 
#'   prior to consensus clustering by \code{cons.funs}. Defaults to 
#'   \code{FALSE}.
#' @param reweigh logical; if \code{TRUE}, algorithms are reweighted based on 
#'   internal validity index performance after trimming. Well-performing 
#'   algorithms are given higher weight prior to consensus clustering by 
#'   \code{cons.funs}. Defaults to \code{FALSE}. Ignored if \code{trim = FALSE}.
#' @param evaluate logical; if \code{TRUE} (default), validity indices are 
#'   returned. Internal validity indices are always computed. If \code{ref.cl} 
#'   is not \code{NULL}, then external validity indices will also be computed.
#' @param plot logical; if \code{TRUE}, \code{\link{graph_all}} is called and 
#'   relevant graphs are outputted. Ignored if \code{evaluate = FALSE}.
#' @param ref.cl reference class; a vector of length equal to the number of 
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
#' ref.cl <- data.frame(initCol = rownames(dat)) %>%
#' tidyr::separate(initCol,
#'                 into = c("patientID", "Class"),
#'                 sep = "_") %>% 
#'   magrittr::use_series(Class) %>% 
#'   factor() %>% 
#'   as.integer()
#' dice.obj <- dice(dat, nk = 4, reps = 5, algorithms = c("hcAEucl", "pamEucl"),
#' cons.funs = c("kmodes", "majority"), ref.cl = ref.cl)
#' str(dice.obj)
dice <- function(data, nk, reps = 10, algorithms = NULL,
                 cons.funs = c("kmodes", "CSPA", "majority", "LCE"),
                 sim.mat = c("cts", "srs", "asrs"),
                 trim = FALSE, reweigh = FALSE, evaluate = TRUE, plot = FALSE,
                 ref.cl = NULL, progress = TRUE) {
  
  # Check that inputs are correct
  assertthat::assert_that(length(dim(data)) == 2)
  if (!is.null(ref.cl)) assertthat::assert_that(is.integer(ref.cl))
  n <- nrow(data)
  ncf <- length(cons.funs)
  
  # Generate Diverse Cluster Ensemble
  E <- consensus_cluster(data = data, nk = nk, reps = reps,
                         algorithms = algorithms, progress = progress)
  
  # Select k
  k <- consensus_evaluate(data = data, E, ref.cl = ref.cl, plot = FALSE)$k
  
  # Evaluate, trim, and reweigh
  if (length(algorithms) > 1 & trim) {
    trim.obj <- consensus_trim(data, E, ref.cl = ref.cl, reweigh = reweigh)
    E <- trim.obj$data.new
  }
  
  # Impute Missing Values using KNN and majority vote
  imp.obj <- impute_missing(E, data)
  Ecomp <- imp.obj$complete
  
  # Consensus functions
  Final <- sapply(cons.funs, function(x) {
    switch(x,
           kmodes = k_modes(Ecomp),
           majority = majority_voting(Ecomp),
           CSPA = CSPA(E, k),
           LCE = LCE(drop(Ecomp), k = k,
                     sim.mat = match.arg(sim.mat))
    )
  })

  # Relabel Final Clustering using reference
  if (is.null(ref.cl)) {
    # Don't relabel if only one consensus function and no reference class
    FinalR <- apply(Final, 2, as.integer)
  } else {
    # Final classes need to be integer for certain functions to work
    FinalR <- apply(Final, 2, function(x) as.integer(
      relabel_class(pred.cl = x, ref.cl = ref.cl)))
  }
  
  # Return evaluation output including consensus function results
  if (evaluate)
    eval.obj <- consensus_evaluate(data, E, cons.cl = FinalR, ref.cl = ref.cl,
                                   plot = plot)
  
  # Add the reference class as the first column if provided
  if (!is.null(ref.cl)) {
    FinalR <- cbind(ref.cl, FinalR)
  }
  rownames(FinalR) <- rownames(data)
  
  return(list(clusters = FinalR, indices = eval.obj))
}
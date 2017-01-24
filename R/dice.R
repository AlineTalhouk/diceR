#' Diverse Clustering Ensemble
#' 
#' Runs consensus clustering across subsamples, algorithms, and number of 
#' clusters (k).
#' 
#' The \code{min.sd} argument is used to filter the feature space for only 
#' highly variable features. Only features with a standard deviation across all 
#' samples greater than \code{min.sd} will be used.
#' 
#' @param data matrix with rows as observations, columns as variables
#' @param nk number of clusters (k) requested; can specify a single integer or a
#'   range of integers to compute multiple k
#' @param reps number of data subsamples to generate. See 
#'   \code{\link{consensus_cluster}} for details.
#' @param algorithms clustering algorithms to be used in the ensemble. Current 
#'   options are "nmf", "hc", "diana", "km", "pam", "ap", "sc", "gmm", "block".
#'   See \code{\link{consensus_cluster}} for details.
#' @param nmf.method specify NMF-based algorithms to run. By default the 
#'   "brunet" and "lee" algorithms are called. See
#'   \code{\link{consensus_cluster}} for details.
#' @param distance a vector of distance functions. Defaults to "euclidean". Can 
#'   use a custom distance function. See \code{\link{consensus_cluster}} for
#'   details.
#' @param cons.funs consensus functions to use. Current options are "kmodes" 
#'   (k-modes), "majority" (majority voting), "CSPA" (Cluster-based Similarity 
#'   Partitioning Algorithm), "LCE" (linkage clustering ensemble)
#' @param sim.mat type of similarity matrix. One of "cts", "srs", "asrs. See 
#'   \code{\link{LCE}} for details.
#' @param min.sd minimum standard deviation threshold. See details.
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
#' dice.obj <- dice(dat, nk = 4, reps = 3, algorithms = "hc", cons.funs =
#' "kmodes", ref.cl = ref.cl)
#' str(dice.obj, max.level = 2)
dice <- function(data, nk, reps = 10, algorithms = NULL,
                 nmf.method = c("brunet", "lee"), distance = "euclidean",
                 cons.funs = c("kmodes", "CSPA", "majority", "LCE"),
                 sim.mat = c("cts", "srs", "asrs"), min.sd = 1,
                 trim = FALSE, reweigh = FALSE, evaluate = TRUE, plot = FALSE,
                 ref.cl = NULL, progress = TRUE) {
  
  # Check that inputs are correct
  assertthat::assert_that(length(dim(data)) == 2)
  n <- nrow(data)
  ncf <- length(cons.funs)
  data <- prepare_data(data, min.sd = min.sd)
  
  # Generate Diverse Cluster Ensemble
  E <- consensus_cluster(data = data, nk = nk, reps = reps,
                         algorithms = algorithms, nmf.method = nmf.method,
                         distance = distance, progress = progress)
  
  # Select k
  k <- consensus_evaluate(data = data, E, ref.cl = ref.cl, plot = FALSE)$k
  
  # Evaluate, trim, and reweigh
  if (length(algorithms) > 1 & trim) {
    trim.obj <- consensus_trim(data, E, ref.cl = ref.cl, reweigh = reweigh)
    E <- trim.obj$data.new
  }
  
  # Impute Missing Values using KNN and majority vote
  imp.obj <- impute_missing(E, data, k)
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
  }) %>%
    apply(2, as.integer)

  # Relabel Final Clustering using reference
  # Don't relabel if only one consensus function and no reference class
  if (!is.null(ref.cl)) {
    FinalR <- Final %>% 
      apply(2, relabel_class, ref.cl = ref.cl)
  } else {
    FinalR <- Final
  }
  
  # Return evaluation output including consensus function results
  if (evaluate)
    eval.obj <- consensus_evaluate(data, E, cons.cl = FinalR, ref.cl = ref.cl,
                                   plot = plot)
  
  # Add the reference class as the first column if provided
  if (!is.null(ref.cl)) {
    FinalR <- cbind(Reference = ref.cl, FinalR)
  }
  rownames(FinalR) <- rownames(data)
  
  return(list(clusters = FinalR, indices = eval.obj))
}

#' Prepare data for consensus clustering
#'
#' Remove variables with low signal and scale before consensus clustering
#'
#' The \code{min.sd} argument is used to filter the feature space for only
#' highly variable features. Only features with a standard deviation across all
#' samples greater than \code{min.sd} will be used.
#'
#' @param data data matrix with rows as samples and columns as variables
#' @param min.sd minimum standard deviation threshold. See details.
#' @return dataset prepared for usage in \code{consensus_cluster}
#' @author Derek Chiu
#' @export
#' @examples
#' set.seed(2)
#' x <- replicate(10, rnorm(100))
#' x.prep <- prepare_data(x)
#' dim(x)
#' dim(x.prep)
prepare_data <- function(data, min.sd = 1) {
  dat.out <- data %>%
    magrittr::extract(apply(., 1, function(x) !any(is.na(x))),
                      apply(., 2, function(x) stats::sd(x, na.rm = TRUE)) >
                        min.sd) %>%
    scale()
  return(dat.out)
}
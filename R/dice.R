#' Diverse Clustering Ensemble
#' 
#' Runs consensus clustering across subsamples, algorithms, and number of 
#' clusters (k).
#' 
#' There are three ways to handle the input data before clustering via argument
#' \code{prep.data}. The default is to use the raw data as-is ("none"). Or, we
#' can enact \code{\link{prepare_data}} on the full dataset ("full"), or the
#' bootstrap sampled datasets ("sampled").
#' 
#' @param data matrix with rows as observations, columns as variables
#' @param nk number of clusters (k) requested; can specify a single integer or a
#'   range of integers to compute multiple k
#' @param reps number of data subsamples to generate. See 
#'   \code{\link{consensus_cluster}} for details.
#' @param algorithms clustering algorithms to be used in the ensemble. Current 
#'   options are "nmf", "hc", "diana", "km", "pam", "ap", "sc", "gmm", "block".
#'   See \code{\link{consensus_cluster}} for details.
#' @param k.method how is k chosen? The default is to use the PAC to choose the 
#'   best k. Specifying an integer as a user-desired k will override the best k 
#'   chosen by PAC. Finally, specifying "all" will produce consensus
#'   results for all k. The "all" method is implicitly performed when the 
#'   number of k used is one.
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
#' @param prep.data Prepare the data on the "full" dataset, the
#'   "sampled" dataset, or "none" (default). See Details.
#' @param min.var minimum variability measure threshold. See
#'   \code{\link{prepare_data}}.
#' @param seed seed used for imputation
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
#' @inheritParams consensus_evaluate
#' @return A list with the following elements
#' \item{E}{raw clustering ensemble object}
#' \item{Eknn}{clustering ensemble object with knn imputation used on \code{E}}
#' \item{Ecomp}{flattened ensemble object with remaining missing entries imputed
#' by majority voting}
#' \item{clusters}{final clustering assignment from the diverse clustering
#' ensemble method}
#' \item{indices}{if \code{evaluate = TRUE}, shows cluster evaluation indices;
#' otherwise \code{NULL}}
#' @author Aline Talhouk, Derek Chiu
#' @export
#' @examples 
#' library(dplyr)
#' data(hgsc)
#' dat <- t(hgsc[, -1])[1:100, 1:50]
#' ref.cl <- data.frame(initCol = rownames(dat)) %>%
#' tidyr::separate(initCol,
#'                 into = c("patientID", "Class"),
#'                 sep = "_") %>%
#'   magrittr::use_series(Class) %>%
#'   factor() %>%
#'   as.integer()
#' dice.obj <- dice(dat, nk = 4, reps = 5, algorithms = "hc", cons.funs =
#' "kmodes", ref.cl = ref.cl, progress = FALSE)
#' str(dice.obj, max.level = 2)
dice <- function(data, nk, reps = 10, algorithms = NULL, k.method = NULL,
                 nmf.method = c("brunet", "lee"), distance = "euclidean",
                 cons.funs = c("kmodes", "CSPA", "majority", "LCE"),
                 sim.mat = c("cts", "srs", "asrs"),
                 prep.data = c("none", "full", "sampled"), min.var = 1,
                 seed = 1,
                 trim = FALSE, reweigh = FALSE, n = 5, 
                 evaluate = TRUE, plot = FALSE,
                 ref.cl = NULL, progress = TRUE) {
  
  # Check that inputs are correct
  assertthat::assert_that(length(dim(data)) == 2)
  prep.data <- match.arg(prep.data)
  
  # Generate Diverse Cluster Ensemble
  E <- consensus_cluster(data = data, nk = nk, reps = reps,
                         algorithms = algorithms, nmf.method = nmf.method,
                         distance = distance, prep.data = prep.data,
                         min.var = min.var, progress = progress)
  # KNN imputation
  Eknn <- apply(E, 2:4, impute_knn, data = data, seed = seed)
  
  # Select k and new (trimmed and reweighed) data
  eval.obj <- consensus_evaluate(data = data, Eknn, ref.cl = ref.cl,
                                 k.method = k.method, trim = trim,
                                 reweigh = reweigh, n = n)
  Eknn <- eval.obj$trim$data.new
  k <- eval.obj$k

  # Impute remaining missing cases
  Ecomp <- impute_missing(Eknn, data, k)
  
  # Consensus functions
  Final <- vapply(cons.funs, function(x) {
    switch(x,
           kmodes = k_modes(Ecomp),
           majority = majority_voting(Ecomp),
           CSPA = CSPA(E, k),
           LCE = LCE(drop(Ecomp), k = k, sim.mat = sim.mat)
    )
  }, double(nrow(Ecomp))) %>%
    apply(2, as.integer)

  # Relabel Final Clustering using reference
  # Don't relabel if only one consensus function and no reference class
  if (is.null(ref.cl)) {
    if (length(cons.funs) == 1) {
      FinalR <- Final
      # If no reference class, > 1 consensus function, use Final[, 1] as ref.cl
    } else {
      FinalR <- apply(Final, 2, relabel_class, ref.cl = Final[, 1])
    }
  } else {
    FinalR <- apply(Final, 2, relabel_class, ref.cl = ref.cl)
  }
  
  # Return evaluation output including consensus function results
  if (evaluate) {
    eval.obj <- consensus_evaluate(data, E, cons.cl = FinalR, ref.cl = ref.cl,
                                   plot = plot)
  } else {
    eval.obj <- NULL
  }

  # Add the reference class as the first column if provided
  if (!is.null(ref.cl)) {
    FinalR <- cbind(Reference = ref.cl, FinalR)
  }
  rownames(FinalR) <- rownames(data)
  
  return(list(E = E, Eknn = Eknn, Ecomp = Ecomp,
              clusters = FinalR, indices = eval.obj))
}

#' Prepare data for consensus clustering
#' 
#' Remove variables with low signal and (optionally) scale before consensus
#' clustering.
#' 
#' @param data data matrix with rows as samples and columns as variables
#' @param scale logical; should the data be centered and scaled?
#' @param type if we use "conventional" measures (default), then the mean and
#'   standard deviation are used for centering and scaling, respectively. If
#'   "robust" measures are specified, the median and median absolute deviation
#'   (MAD) are used.
#' @param min.var minimum variability measure threshold used to filter the
#'   feature space for only highly variable features. Only features with a
#'   minimum variability measure across all samples greater than \code{min.var}
#'   will be used. If \code{type = "conventional"}, the standard deviation is
#'   the measure used, and if \code{type = "robust"}, the MAD is the measure
#'   used.
#' @return dataset prepared for usage in \code{consensus_cluster}
#' @author Derek Chiu
#' @export
#' @examples
#' set.seed(2)
#' x <- replicate(10, rnorm(100))
#' x.prep <- prepare_data(x)
#' dim(x)
#' dim(x.prep)
prepare_data <- function(data, scale = TRUE,
                         type = c("conventional", "robust"),
                         min.var = 1) {
  type <- match.arg(type)
  var.fun <- switch(type, conventional = stats::sd, robust = stats::mad)
  dat <- data %>% 
    magrittr::extract(stats::complete.cases(.),
                      apply(., 2, var.fun, na.rm = TRUE) > min.var)
  if (scale) {
    dat <- switch(type,
                   conventional = scale(dat),
                   robust = quantable::robustscale(dat))
  }
  return(dat)
}
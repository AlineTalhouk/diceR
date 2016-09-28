#' Run a series of ConsensusClusterPlus algorithms
#'
#' Runs consensus clustering for multiple algorithms and saves the output
#'
#' @param dat data matrix; columns are samples and rows are genes/features
#' @param k maximum number of clusters to compute consensus results
#' @param reps number of subsamples
#' @param pItem proportion of features to use in each subsample
#' @param save logical; if \code{TRUE}, the returned object will be saved at
#'   each iteration as well as at the end.
#' @param file.name file name of the written object
#' @param seed random seed to maintain reproducible results
#' @param min.sd minimum standard deviation across each gene
#' @param ... additional arguments to \code{ConsensusClusterPlus}
#' @return A list of outputs from \code{ConsensusClusterPlus}; each element
#' is for a different algorithm.
#' @author Derek Chiu
#' @import ConsensusClusterPlus
#' @export
#' @examples
#' set.seed(23)
#' x <- ConClustPlus(matrix(rnorm(100), nrow = 10), k = 3, reps = 10, pItem = 0.9)
ConClustPlus <- function(dat, k = 3, reps = 1000, pItem = 0.8, save = FALSE,
                         file.name = "results_CCP", seed = 123, min.sd = 1,
                         ...) {
  .Deprecated("ConClust")
  dat <- prepare_data(dat, min.sd = min.sd)
  algs <- setNames(c("hc", "diana_hook", "kmdist", "kmdist", "pam", "pam"),
                   c("hcAEucl", "hcDianaEucl", "kmEucl", "kmSpear", "pamEucl",
                     "pamSpear"))
  dists <- c("euclidean", "euclidean", "euclidean", "spearman", "euclidean",
             "spearman")
  results <- mapply(function(a, d)
    ConsensusClusterPlus(dat, maxK = k, reps = reps, pItem = pItem,
                         clusterAlg = a, distance = d,
                         seed = seed, ...),
    a = algs, d = dists, SIMPLIFY = FALSE)
  if (save)
    readr::write_rds(results, paste0(file.name, ".rds"), compress = "xz")
  return(results)
}

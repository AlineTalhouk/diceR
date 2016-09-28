#' Combine and compare consensus results
#'
#' Combine results from ConClust and ConClustPlus and output either the
#' consensus matrices or consensus classes for all algorithms from both objects.
#' Compare algorithms on validation indices PAC and CHI and weigh algorithms
#' based on these two measures.
#'
#' @param ... any number of objects outputted from
#'   \code{\link{consensus_summary}}
#' @param res.CCP an object outputted from \code{\link{ConClustPlus}}
#' @param k desired number of clusters
#' @param element either "matrix" or "class" to extract the consensus matrix or
#'   consensus class, respectively.
#' @param alg.names optional. Supply a vector of names for the algorithms.
#' @return \code{consensus_combine} returns either a list of all consensus
#'   matrices or a data frame showing all the consensus classes
#' @family consensus functions
#' @author Derek Chiu
#' @export
#' @examples
#' # Consensus clustering for multiple algorithms
#' set.seed(911)
#' x <- matrix(rnorm(1000), nrow = 10)
#' CC1 <- ConClust(x, k = 4, reps = 10, method = "apEucl", save = FALSE)
#' CC2 <- ConClust(x, k = 4, reps = 10, method = "gmmBIC", save = FALSE)
#' CCP <- ConClustPlus(x, k = 4, reps = 10, save = FALSE)
#'
#' # Get summary for ConClust
#' CC1.summ <- consensus_summary(CC1, k = 4)
#' CC2.summ <- consensus_summary(CC2, k = 4)
#'
#' # Combine with CCP and return either matrices or classes
#' y1 <- consensus_combine(CC1.summ, CC2.summ, res.CCP = CCP, k = 4,
#' element = "matrix")
#' str(y1)
#' y2 <- consensus_combine(CC1.summ, CC2.summ, res.CCP = CCP, k = 4,
#' element = "class")
#' str(y2)
#'
#' # Compare algorithms on PAC and CHI
#' z <- consensus_compare(x, cl.mat = y2, cons.mat = y1)
#'
#' # Weigh algorithms
#' consensus_weigh(z)
consensus_combine <- function(..., res.CCP, k, element = c("matrix", "class"),
                              alg.names = NULL) {
  obj <- unlist(list(...), recursive = FALSE)
  switch(match.arg(element),
         matrix = {
           out.CC <- lapply(obj, "[[", "consensus_matrix")
           out.CCP <- lapply(lapply(res.CCP, "[[", k), "[[", "consensusMatrix")
           out.CCP <- lapply(out.CCP, function(x) {
             dimnames(x) <- dimnames(out.CC[[1]])
             return(x)
           })
           out <- unlist(list(out.CCP, out.CC), recursive = FALSE)
           if (!is.null(alg.names))
             names(out) <- alg.names
         },
         class = {
           out.CC <- apply(sapply(obj, "[[", "consensus_class"), c(1, 2),
                           as.integer)
           out.CCP <- sapply(lapply(res.CCP, "[[", k), "[[", "consensusClass")
           out <- as.matrix(data.frame(out.CCP, out.CC))
           if (!is.null(alg.names))
             colnames(out) <- alg.names
         })
  return(out)
}

#' Combine, compare, and weigh algorithms
#'
#' \code{consensus_combine} combines results for multiple objects from
#' \code{consensus_summary(ConClust())} and outputs either the consensus
#' matrices or consensus classes for all algorithms. \code{consensus_compare}
#' compares algorithms on validation indices \code{\link{PAC}} and CHI.
#' \code{consensus_weigh} weighs clustering algorithms based on these two
#' indices.
#'
#' \code{consensus_combine} is useful for generating summaries because the
#' results have been combined into a single object. For example, if
#' \code{element = "class"}, then the resulting object can be used to create a
#' consensus matrix across algorithms, which can be visualized as a heatmap.
#'
#' @param ... any number of objects outputted from
#'   \code{\link{consensus_summary}}
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
#'
#' # Get summary for ConClust
#' CC1.summ <- consensus_summary(CC1, k = 4)
#' CC2.summ <- consensus_summary(CC2, k = 4)
#'
#' # Combine and return either matrices or classes
#' y1 <- consensus_combine(CC1.summ, CC2.summ, element = "matrix")
#' str(y1)
#' y2 <- consensus_combine(CC1.summ, CC2.summ, element = "class")
#' str(y2)
#'
#' # Compare algorithms on PAC and CHI
#' z <- consensus_compare(x, cl.mat = y2, cons.mat = y1)
#'
#' # Weigh algorithms
#' consensus_weigh(z)
consensus_combine <- function(..., element = c("matrix", "class"),
                              alg.names = NULL) {
  obj <- unlist(list(...), recursive = FALSE)
  switch(match.arg(element),
         matrix = {
           out <- lapply(obj, "[[", "consensus_matrix")
           out <- unlist(list(out), recursive = FALSE)
           if (!is.null(alg.names))
             names(out) <- alg.names
         },
         class = {
           out <- apply(sapply(obj, "[[", "consensus_class"), c(1, 2),
                        as.integer)
           out <- as.matrix(data.frame(out))
           if (!is.null(alg.names))
             colnames(out) <- alg.names
         })
  return(out)
}

#' Combine, evaluate, trim, and reweigh algorithms
#'
#' \code{consensus_combine} combines results for multiple objects from
#' \code{ConClust()} and outputs either the consensus
#' matrices or consensus classes for all algorithms. \code{consensus_evaluate}
#' evaluates algorithms on internal/external validation indices.
#' \code{consensus_weigh} weighs clustering algorithms based on these two
#' indices. \code{consensus_trim} removes algorithms that rank low on internal
#' indices before using in ensemble clustering methods.
#'
#' \code{consensus_combine} is useful for generating summaries because the
#' results have been combined into a single object. For example, if
#' \code{element = "class"}, then the resulting object can be used to create a
#' consensus matrix across algorithms, which can be visualized as a heatmap.
#' 
#' \code{consensus_evaluate} always shows internal indices. If \code{ref.cl} is
#' not \code{NULL}, external indices are shown in addition to internal indices.
#' Relevant graphical displays are also outputted.
#' 
#' \code{consensus_trim} ranks algorithms by internal indices used in
#' \code{consensus_evaluate}. The sum of the ranks for each algorithm is used as
#' the measure of comparison. This also means the magnitude of the internal
#' indices is not taken into account.
#'
#' @param ... any number of objects outputted from
#'   \code{\link{ConClust}}
#' @param k number of clusters requested. Default is \code{NULL}, which returns
#'   all k computed from \code{...}
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
#' x <- matrix(rnorm(500), ncol = 10)
#' CC1 <- ConClust(x, nc = 3:4, reps = 10, method = "apEucl")
#' CC2 <- ConClust(x, nc = 3:4, reps = 10, method = "gmmBIC")
#' 
#' # Combine and return either matrices or classes
#' y1 <- consensus_combine(CC1, CC2, element = "matrix")
#' str(y1)
#' y2 <- consensus_combine(CC1, CC2, element = "class")
#' str(y2)
#' 
#' # Evaluate algorithms on internal and external indices and make plots
#' set.seed(1)
#' ref.cl <- sample(1:4, 50, replace = TRUE)
#' z <- consensus_evaluate(x, CC1, CC2, ref.cl = ref.cl, plot = FALSE)
#' 
#' # Trim algorithms: remove those that rank low on internal indices
#' CC3 <- consensus_trim(x, CC1, CC2, ref.cl = ref.cl, quantile = 0.8)
#' str(CC3)
consensus_combine <- function(..., element = c("matrix", "class"),
                              alg.names = NULL) {
  cs <- abind::abind(list(...), along = 3)
  obj <- consensus_summary(cs)
  switch(match.arg(element),
         matrix = {
           out <- lapply(obj, purrr::transpose) %>% 
             lapply("[[", "consensus_matrix")
           if (!is.null(alg.names))
             out <- lapply(out, function(x) x %>% set_names(alg.names))
         },
         class = {
           out <- lapply(obj, purrr::transpose) %>% 
             lapply("[[", "consensus_class") %>% 
             lapply(as.data.frame) %>% 
             lapply(function(x) apply(x, 1:2, as.integer))
           if (!is.null(alg.names))
             out <- lapply(out, function(x) x %>% set_colnames(alg.names))
         })
  return(out)
}
#' Combine, evaluate, trim, and reweigh algorithms
#'
#' \code{consensus_combine} combines results for multiple objects from
#' \code{consensus_cluster()} and outputs either the consensus
#' matrices or consensus classes for all algorithms. \code{consensus_evaluate}
#' evaluates algorithms on internal/external validation indices.
#' \code{consensus_weigh} weighs clustering algorithms based on these two
#' indices. \code{consensus_trim} removes algorithms that rank low on internal
#' indices before using in consensus functions.
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
#'   \code{\link{consensus_cluster}}
#' @param element either "matrix" or "class" to extract the consensus matrix or
#'   consensus class, respectively.
#' @return \code{consensus_combine} returns either a list of all consensus
#'   matrices or a data frame showing all the consensus classes
#' @author Derek Chiu
#' @export
#' @examples
#' # Consensus clustering for multiple algorithms
#' set.seed(911)
#' x <- matrix(rnorm(500), ncol = 10)
#' CC1 <- consensus_cluster(x, nk = 3:4, reps = 10, algorithms = "ap",
#' progress = FALSE)
#' CC2 <- consensus_cluster(x, nk = 3:4, reps = 10, algorithms = "gmm",
#' progress = FALSE)
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
#' str(CC3, max.level = 2)
consensus_combine <- function(..., element = c("matrix", "class")) {
  # Combine ensemble arrays and reorganize into matrices and classes
  cs <- abind::abind(list(...), along = 3)
  # Return a list of summaries for each algorithm
  obj <- consensus_summary(cs)
  switch(match.arg(element),
         matrix = {
           # Transpose list levels and extract matrices
           out <- lapply(obj, purrr::transpose) %>% 
             lapply("[[", "consensus_matrix")
         },
         class = {
           # Transpose list levels and extract classes, coercing to integer
           out <- lapply(obj, purrr::transpose) %>% 
             lapply("[[", "consensus_class") %>% 
             lapply(as.data.frame) %>% 
             lapply(function(x) apply(x, 1:2, as.integer))
         })
  return(out)
}

#' Given an object from \code{\link{consensus_cluster}}, returns a list of
#' consensus matrices and consensus classes for each clustering algorithm.
#' @noRd
consensus_summary <- function(E) {
  con.mats <- plyr::alply(E, 3:4, consensus_matrix, .dims = TRUE) %>% 
    utils::relist(stats::setNames(replicate(
      dim(E)[4],
      list(structure(1:dim(E)[3], names = dimnames(E)[[3]]))),
      dimnames(E)[[4]]))
  con.cls <- mapply(function(cm, k) lapply(cm, function(x) hc(stats::dist(x),
                                                              k = k)),
                    cm = con.mats, k = as.numeric(names(con.mats)),
                    SIMPLIFY = FALSE)
  out <- list(consensus_matrix = con.mats, consensus_class = con.cls) %>% 
    purrr::transpose() %>% 
    lapply(purrr::transpose)
  return(out)
}
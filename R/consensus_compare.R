#' @param data data matrix with samples as rows, genes as columns
#' @param cl.mat matrix of cluster assignments. Each row is an assignment for a
#'   different algorithm. Use \code{element = "class"} in
#'   \code{consensus_combine}.
#' @param cons.mat A list of consensus matrices, one for each algorithm. Use
#'   \code{element = "matrix"} in \code{consensus_combine}.
#' @inheritParams consensus_combine
#' @return \code{consensus_compare} returns a data frame of the indices PAC and
#'   CHI in each column for each algorithm.
#' @rdname consensus_combine
#' @export
consensus_compare <- function(data, cl.mat, cons.mat, alg.names = NULL) {
  if (!is.null(alg.names)) {
    an <- alg.names
  } else {
    an <- names(cons.mat)
  }
  ind <- data.frame(Algorithms = an,
                    PAC = sapply(cons.mat, PAC, lower = 0.05, upper = 0.95),
                    CHI = apply(cl.mat, 2, clusterSim::index.G1, x = t(data)))
  return(ind)
}

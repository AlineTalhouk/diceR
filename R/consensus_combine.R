#' Combine algorithms
#'
#' Combines results for multiple objects from \code{consensus_cluster()} and
#' outputs either the consensus matrices or consensus classes for all
#' algorithms.
#'
#' This function is useful for collecting summaries because the original results
#' from \code{consensus_cluster} were combined to a single object. For example,
#' setting \code{element = "class"} returns a matrix of consensus cluster
#' assignments, which can be visualized as a consensus matrix heatmap.
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
consensus_combine <- function(..., element = c("matrix", "class")) {
  # Combine ensemble arrays and reorganize into matrices and classes
  cs <- abind::abind(list(...), along = 3)
  # Return a list of summaries for each algorithm
  obj <- consensus_summary(cs)
  switch(match.arg(element),
         matrix = {
           # Transpose list levels and extract matrices
           out <- purrr::map(obj, purrr::transpose) %>%
             purrr::map("consensus_matrix")
         },
         class = {
           # Transpose list levels and extract classes, coercing to integer
           out <- purrr::map(obj, purrr::transpose) %>%
             purrr::map("consensus_class") %>%
             purrr::map(as.data.frame) %>%
             purrr::map(~ apply(.x, 1:2, as.integer))
         })
  out
}

#' Given an object from \code{\link{consensus_cluster}}, returns a list of
#' consensus matrices and consensus classes for each clustering algorithm.
#' @noRd
consensus_summary <- function(E) {
  con.mats <- E %>%
    purrr::array_tree(c(4, 3)) %>%
    purrr::at_depth(2, consensus_matrix) %>%
    purrr::map(magrittr::set_names, dimnames(E)[[3]]) %>%
    magrittr::set_names(dimnames(E)[[4]])
  con.cls <- purrr::map2(con.mats, as.numeric(names(con.mats)),
                         ~ purrr::map(.x, function(z)
                           hc(stats::dist(z), k = .y)))
  out <- list(consensus_matrix = con.mats, consensus_class = con.cls) %>%
    purrr::transpose() %>%
    purrr::map(purrr::transpose)
  out
}

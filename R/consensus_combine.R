#' Combine algorithms
#'
#' Combines results for multiple objects from `consensus_cluster()` and outputs
#' either the consensus matrices or consensus classes for all algorithms.
#'
#' This function is useful for collecting summaries because the original results
#' from `consensus_cluster` were combined to a single object. For example,
#' setting `element = "class"` returns a matrix of consensus cluster
#' assignments, which can be visualized as a consensus matrix heatmap.
#'
#' @param ... any number of objects outputted from [consensus_cluster()]
#' @param element either "matrix" or "class" to extract the consensus matrix or
#'   consensus class, respectively.
#' @return `consensus_combine` returns either a list of all consensus matrices
#'   or a data frame showing all the consensus classes
#' @author Derek Chiu
#' @export
#' @examplesIf rlang::is_installed("apcluster")
#' # Consensus clustering for multiple algorithms
#' set.seed(911)
#' x <- matrix(rnorm(500), ncol = 10)
#' CC1 <- consensus_cluster(x, nk = 3:4, reps = 10, algorithms = "ap",
#' progress = FALSE)
#' CC2 <- consensus_cluster(x, nk = 3:4, reps = 10, algorithms = "km",
#' progress = FALSE)
#'
#' # Combine and return either matrices or classes
#' y1 <- consensus_combine(CC1, CC2, element = "matrix")
#' str(y1)
#' y2 <- consensus_combine(CC1, CC2, element = "class")
#' str(y2)
consensus_combine <- function(..., element = c("matrix", "class")) {
  cs <- abind::abind(list(...), along = 3) %>% # Bind ensemble arrays on algs
    consensus_summary() # Reorganize into matrices and classes
  switch(
    match.arg(element),
    matrix = purrr::map(cs, "con.mats"),
    class = purrr::map(cs, "con.cls") %>%
      purrr::map(~ do.call(cbind, .)) # Combine classes into list of matrices
  )
}

#' Given an object from [consensus_cluster()], returns a list of consensus
#' matrices and consensus classes for each clustering algorithm.
#' @noRd
consensus_summary <- function(E) {
  con.mats <- E %>%
    purrr::array_tree(c(4, 3)) %>%
    purrr::modify_depth(2, consensus_matrix) %>%
    purrr::map(magrittr::set_names, dimnames(E)[[3]]) %>%
    magrittr::set_names(dimnames(E)[[4]])
  con.cls <- con.mats %>%
    purrr::imap(~ purrr::map(.x, function(z) hc(stats::dist(z), k = .y)))
  dplyr::lst(con.mats, con.cls) %>% purrr::transpose() # transpose lists
}

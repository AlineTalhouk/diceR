#' @param data data matrix with rows as samples and columns as variables
#' @param ref.cl reference class
#' @param cons.cl matrix of cluster assignments from consensus algorithms such
#'   as \code{kmodes} and \code{majority_voting}
#' @param plot logical; if \code{TRUE}, \code{graph_all} is called
#' @inheritParams consensus_combine
#' @return \code{consensus_evaluate} returns a data frame of the indices in each
#'   column for each algorithm.
#' @rdname consensus_combine
#' @export
consensus_evaluate <- function(data, k = NULL, ..., cons.cl = NULL,
                               ref.cl = NULL, plot = TRUE) {
  x <- as.matrix(data)
  cc.obj <- abind::abind(list(...), along = 3)
  # cl.mat <- consensus_combine(cc.obj, k = k, element = "class")
  cl.mat <- consensus_combine(cc.obj, element = "class")
  an <- dimnames(cc.obj)[3][[1]]
  if (!is.null(cons.cl)) {
    assertthat::assert_that(is.matrix(cons.cl))
    # cl.mat <- cbind(cl.mat, cons.cl)
    cl.mat <- lapply(cl.mat, cbind, cons.cl)
    an <- c(an, colnames(cons.cl))
  }
  # cl.mat <- apply(cl.mat, 1:2, as.integer)
  
  ind.int <- lapply(cl.mat, function(m) {
    data.frame(
      Algorithms = an,
      plyr::aaply(m, 2, function(cl)
        clusterCrit::intCriteria(
          traj = x, part = cl,
          crit = c("C_index", "Calinski_Harabasz",
                   "Davies_Bouldin", "Dunn", "McClain_Rao",
                   "PBM", "SD_Dis", "Ray_Turi", "Tau",
                   "Gamma", "G_plus")) %>%
          unlist()),
      Compactness = apply(m, 2, iv_compactness, data = x),
      Connectivity = apply(m, 2, function(cl)
        clValid::connectivity(Data = x, clusters = cl))) %>%
      mutate_all(funs(structure(., names = an)))
  })
  # ind.int <- data.frame(
  #   Algorithms = an,
  #   plyr::aaply(cl.mat, 2, function(cl)
  #     clusterCrit::intCriteria(
  #       traj = x, part = cl,
  #       crit = c("C_index", "Calinski_Harabasz",
  #                "Davies_Bouldin", "Dunn", "McClain_Rao",
  #                "PBM", "SD_Dis", "Ray_Turi", "Tau",
  #                "Gamma", "G_plus")) %>%
  #       unlist()),
  #   Compactness = apply(cl.mat, 2, iv_compactness, data = x),
  #   Connectivity = apply(cl.mat, 2, function(cl)
  #     clValid::connectivity(Data = x, clusters = cl))) %>%
  #   mutate_all(funs(structure(., names = an)))
  if (plot) {
    graph_all(cc.obj)
  }
  if (!is.null(ref.cl)) {
    cl.mat.ext <- cl.mat %>% 
      extract2(match(n_distinct(ref.cl), names(.)))
    ind.ext <- data.frame(
      Algorithms = an,
      plyr::aaply(cl.mat.ext, 2, function(cl)
        clusterCrit::extCriteria(
          part1 = cl, part2 = ref.cl,
          crit = c("Hubert", "Jaccard", "McNemar",
                   "Precision", "Rand", "Recall")) %>%
          unlist()),
      NMI = apply(cl.mat.ext, 2, ev_nmi, ref.lab = ref.cl)) %>%
      cbind(t(apply(cl.mat.ext, 2, ev_confmat, ref.lab = ref.cl))) %>%
      mutate_all(funs(structure(., names = an)))
    return(list(internal = ind.int, external = ind.ext))
  } else {
    return(list(internal = ind.int))
  }
}

#' Compactness Measure
#' 
#' Compute the compactness validity index for a clustering result.
#' 
#' This index is agnostic to any reference clustering results, calculating
#' cluster performance on the basis of compactness and separability.
#' 
#' @param data a dataset with rows as observations, columns as variables
#' @param labels a vector of cluster labels from a clustering result
#' @references MATLAB function \code{valid_compactness} by Simon Garrett in
#'   LinkCluE
#' @return the compactness score
#' @export
#' 
#' @examples
#' set.seed(1)
#' E <- matrix(rep(sample(1:4, 1000, replace = TRUE)), nrow = 100, byrow =
#'               FALSE)
#' set.seed(1)
#' dat <- as.data.frame(matrix(runif(1000, -10, 10), nrow = 100, byrow = FALSE))
#' iv_compactness(dat, E[, 1])
iv_compactness <- function(data, labels) {
  assertthat::assert_that(is.data.frame(data) || is.matrix(data),
                          length(labels) == nrow(data))
  n <- length(labels)
  C <- sort(unique(labels))
  k <- length(C)
  cp <- 0
  for (i in 1:k) {
    ind <- which(labels == C[i])
    nk <- length(ind)
    if (nk <= 1) {
      cp <- cp + 0
    } else{
      sum_d <- 0
      sum_d <- sum(stats::dist(data[ind, ], method = "euclidean"))
      cp <- cp + (nk * (sum_d / (nk * (nk - 1) / 2)))
    }
  }
  return(cp / n)
}
#' @param data data matrix with rows as samples and columns as variables
#' @param ref.cl reference class
#' @param cons.cl matrix of cluster assignments from consensus functions such
#'   as \code{kmodes} and \code{majority_voting}
#' @param plot logical; if \code{TRUE}, \code{graph_all} is called
#' @inheritParams consensus_combine
#' @return \code{consensus_evaluate} returns a list with the following elements
#'   \item{k}{if \code{ref.cl} is not NULL, this is the number of distinct
#'   classes in the reference; otherwise the chosen \code{k} is determined by
#'   the one giving the largest mean PAC across algorithms}
#'   \item{pac}{a data frame showing the PAC for each combination of algorithm
#'   and cluster size}
#'   \item{internal}{a list of data frames for all k showing internal evaluation
#'   indices}
#'   \item{external}{a data frame showing external evaluation indices for
#'   \code{k}}
#' @rdname consensus_combine
#' @export
consensus_evaluate <- function(data, ..., cons.cl = NULL, ref.cl = NULL,
                               plot = FALSE) {
  # Assertions
  if (!is.null(ref.cl))
    assertthat::assert_that(is.integer(ref.cl), nrow(data) == length(ref.cl))
  
  # Extract classes and matrices separately
  x <- as.matrix(data)
  cc.obj <- abind::abind(list(...), along = 3)
  cl.mat <- consensus_combine(cc.obj, element = "class")
  cons.mat <- consensus_combine(cc.obj, element = "matrix")
  
  # Calculate PAC
  pac <- lapply(cons.mat, lapply, PAC) %>%
    plyr::ldply(unlist, .id = "k")
  
  # If reference given, k is number of distinct classes
  if (!is.null(ref.cl)) {
    k <- n_distinct(ref.cl)
    # Otherwise k is the maximum average PAC across algorithms
  } else {  
    idx.k <- apply(pac[, -1, drop = FALSE], 2, which.max) %>% 
      table() %>% 
      which.max() %>% 
      names() %>% 
      as.numeric()
    k <- as.integer(as.character(pac$k[idx.k]))
  }
  
  # If matrix of cluster assignments from cons.funs given, cbind to cl.mat
  an <- dimnames(cc.obj)[3][[1]]
  if (!is.null(cons.cl)) {
    assertthat::assert_that(is.matrix(cons.cl))
    cl.mat <- lapply(cl.mat, cbind, cons.cl)
    an <- c(an, colnames(cons.cl))
  }
  
  # Internal indices
  ind.int <- lapply(cl.mat, function(m) {
    data.frame(
      Algorithms = an,
      apply(m, 2, function(cl)
        clusterCrit::intCriteria(
          traj = x, part = cl,
          crit = c("C_index", "Calinski_Harabasz",
                   "Davies_Bouldin", "Dunn", "McClain_Rao",
                   "PBM", "SD_Dis", "Ray_Turi", "Tau",
                   "Gamma", "G_plus")) %>%
          unlist()) %>% t(),
      Compactness = apply(m, 2, compactness, data = x),
      Connectivity = apply(m, 2, function(cl)
        clValid::connectivity(Data = x, clusters = cl))) %>%
      mutate_all(funs(structure(., names = an)))
  })
  
  # Graph all plotting functions
  if (plot) {
    graph_all(cc.obj)
  }
  
  # Only calculate external indices if a reference is given
  if (!is.null(ref.cl)) {
    cl.mat.ext <- cl.mat %>% 
      extract2(match(k, names(.)))
    ind.ext <- data.frame(
      Algorithms = an,
      apply(cl.mat.ext, 2, function(cl)
        clusterCrit::extCriteria(
          part1 = cl, part2 = ref.cl,
          crit = c("Hubert", "Jaccard", "McNemar",
                   "Precision", "Rand", "Recall")) %>%
          unlist()) %>% t(),
      NMI = apply(cl.mat.ext, 2, ev_nmi, ref.lab = ref.cl)) %>%
      cbind(t(apply(cl.mat.ext, 2, ev_confmat, ref.lab = ref.cl))) %>%
      mutate_all(funs(structure(., names = an)))
    return(list(k = k, pac = pac, internal = ind.int, external = ind.ext))
  } else {
    return(list(k = k, pac = pac, internal = ind.int))
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
#' compactness(dat, E[, 1])
compactness <- function(data, labels) {
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
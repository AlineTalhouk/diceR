#' @param data data matrix with rows as samples and columns as variables
#' @param cons.cl matrix of cluster assignments from consensus functions such
#'   as \code{kmodes} and \code{majority_voting}
#' @param ref.cl reference class
#' @param k.method determines the method to choose k when no reference class is 
#'   given. When \code{ref.cl} is not \code{NULL}, k is the number of distinct
#'   classes of \code{ref.cl}. Otherwise the input from \code{k.method} chooses
#'   k. The default is to use the PAC to choose the best k(s). Specifying an 
#'   integer as a user-desired k will override the best k chosen by PAC. 
#'   Finally, specifying "all" will produce consensus results for all k. The 
#'   "all" method is implicitly performed when there is only one k used.
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
                               k.method = NULL, plot = FALSE) {
  # Assertions
  if (!is.null(ref.cl))
    assertthat::assert_that(is.integer(ref.cl), nrow(data) == length(ref.cl))
  
  # Extract classes and matrices separately
  x <- as.matrix(data)
  cc.obj <- abind::abind(list(...), along = 3)
  cl.mat <- consensus_combine(cc.obj, element = "class")
  cons.mat <- consensus_combine(cc.obj, element = "matrix")
  
  # Calculate PAC
  pac <- cons.mat %>% 
    purrr::at_depth(2, PAC) %>% 
    purrr::map_df(data.frame, .id = "k")
  
  # If reference given, k is number of distinct classes
  if (!is.null(ref.cl)) {
    k <- n_distinct(ref.cl)
    # Otherwise k is chosen using the following methods
  } else if (is.null(k.method)) {
    k <- pac %>% 
      use_series("k") %>% 
      extract(apply(pac[, -1, drop = FALSE], 2, which.min) %>% 
                table() %>% 
                extract(is_in(., max(.))) %>% 
                names() %>% 
                as.numeric()) %>% 
      as.integer()
  } else if (length(cons.mat) == 1 || k.method == "all") {
    k <- as.integer(dimnames(cc.obj)[[4]])
  } else if (length(k.method) == 1 & is.numeric(k.method)) {
    k <- k.method
  } else {
    stop("Invalid input. Check documentation for possible options.")
  }
  
  # If matrix of cluster assignments from cons.funs given, cbind to cl.mat
  an <- dimnames(cc.obj)[3][[1]]
  if (!is.null(cons.cl)) {
    assertthat::assert_that(is.matrix(cons.cl))
    cl.mat <- purrr::map(cl.mat, cbind, cons.cl)
    an <- c(an, colnames(cons.cl))
  }
  
  # Internal indices
  ind.int <- purrr::map(cl.mat, ~ data.frame(
    Algorithms = an,
    apply(.x, 2, function(cl)
      clusterCrit::intCriteria(
        traj = x, part = cl,
        crit = c("C_index", "Calinski_Harabasz", "Davies_Bouldin", "Dunn",
                 "McClain_Rao", "PBM", "SD_Dis", "Ray_Turi", "Tau", "Gamma",
                 "G_plus", "Silhouette", "S_Dbw")) %>%
        unlist()) %>% t(),
    Compactness = apply(.x, 2, compactness, data = x),
    Connectivity = apply(.x, 2, function(cl)
      clValid::connectivity(Data = x, clusters = cl))) %>%
      mutate_all(funs(structure(., names = an))))
  
  # Graph all plotting functions
  if (plot) {
    graph_all(cc.obj)
  }
  
  # Only calculate external indices if a reference is given
  if (!is.null(ref.cl)) {
    cl.mat.ext <- cl.mat %>% 
      extract(match(k, names(.)))
    ind.ext <- purrr::map(cl.mat.ext, ~ data.frame(
      Algorithms = an,
      apply(.x, 2, function(cl)
        clusterCrit::extCriteria(
          part1 = cl, part2 = ref.cl,
          crit = c("Hubert", "Jaccard", "McNemar", "Precision", "Rand",
                   "Recall")) %>%
          unlist()) %>% t(),
      NMI = apply(.x, 2, ev_nmi, ref.lab = ref.cl)) %>%
        cbind(t(apply(.x, 2, ev_confmat, ref.lab = ref.cl))) %>%
        mutate_all(funs(structure(., names = an))))
  } else {
    ind.ext <- NULL
  }
  return(list(k = k, pac = pac, internal = ind.int, external = ind.ext))
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
#' @return the compactness score
#' @author Derek Chiu
#' @references MATLAB function \code{valid_compactness} by Simon Garrett in
#'   LinkCluE
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
  C <- sort(unique(labels))
  cp <- 0
  for (i in seq_along(C)) {
    ind <- which(labels == C[i])
    nk <- length(ind)
    if (nk > 1) {
      sum_d <- sum(stats::dist(data[ind, ], method = "euclidean"))
      cp <- cp + (nk * (sum_d / (nk * (nk - 1) / 2)))
    }
  }
  return(cp / length(labels))
}
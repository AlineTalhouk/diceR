#' Evaluate, trim, and reweigh algorithms
#'
#' Evaluates algorithms on internal/external validation indices. Poor performing
#' algorithms can be trimmed from the ensemble. The remaining algorithms can be
#' given weights before use in consensus functions.
#'
#' This function always returns internal indices. If \code{ref.cl} is not
#' \code{NULL}, external indices are additionally shown. Relevant graphical
#' displays are also outputted. Algorithms are ranked across internal indices
#' using Rank Aggregation. Only the top \code{n} algorithms are kept, the rest
#' are trimmed.
#'
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
#' @param trim logical; if \code{TRUE}, algorithms that score low on internal
#'   indices will be trimmed out
#' @param reweigh logical; if \code{TRUE}, after trimming out poor performing
#'   algorithms, each algorithm is reweighed depending on its internal indices.
#' @param n an integer specifying the top \code{n} algorithms to keep after
#'   trimming off the poor performing ones using Rank Aggregation. If the total
#'   number of algorithms is less than \code{n} no trimming is done.
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
#'   \item{trim}{A list with 4 elements}
#'   \itemize{
#'     \item{\code{alg.keep} }{algorithms kept}
#'     \item{\code{alg.remove} }{algorithms removed}
#'     \item{\code{rank.agg} }{a matrix of ranked algorithms for every internal
#'     evaluation index}
#'     \item{\code{top.list} }{final order of ranked algorithms}
#'     \item{\code{data.new} }{A new version of a \code{consensus_cluster} data
#'     object}
#'   }
#' @export
#' @examples
#' # Consensus clustering for multiple algorithms
#' set.seed(911)
#' x <- matrix(rnorm(500), ncol = 10)
#' CC <- consensus_cluster(x, nk = 3:4, reps = 10, algorithms = c("ap", "gmm"),
#' progress = FALSE)
#'
#' # Evaluate algorithms on internal/external indices and trim algorithms:
#' # remove those ranking low on internal indices
#' set.seed(1)
#' ref.cl <- sample(1:4, 50, replace = TRUE)
#' z <- consensus_evaluate(x, CC, ref.cl = ref.cl, n = 1, trim = TRUE)
#' str(z, max.level = 2)
consensus_evaluate <- function(data, ..., cons.cl = NULL, ref.cl = NULL,
                               k.method = NULL, plot = FALSE, trim = FALSE,
                               reweigh = FALSE, n = 5) {
  # Assertions
  if (!is.null(ref.cl))
    assertthat::assert_that(is.integer(ref.cl), nrow(data) == length(ref.cl))

  # Extract classes and matrices separately
  E <- abind::abind(list(...), along = 3)
  cl.mat <- consensus_combine(E, element = "class")
  cons.mat <- consensus_combine(E, element = "matrix")

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
      magrittr::extract(apply(pac[, -1, drop = FALSE], 2, which.min) %>%
                          table() %>%
                          magrittr::extract(is_in(., max(.))) %>%
                          names() %>%
                          as.numeric()) %>%
      as.integer() %>%
      min()
  } else if (length(cons.mat) == 1 || k.method == "all") {
    k <- as.integer(dimnames(E)[[4]])
  } else if (length(k.method) == 1 & is.numeric(k.method)) {
    k <- k.method
  } else {
    stop("Invalid input. Check documentation for possible options.")
  }

  # If matrix of cluster assignments from cons.funs given, cbind to cl.mat
  an <- dimnames(E)[[3]]
  if (!is.null(cons.cl)) {
    assertthat::assert_that(is.matrix(cons.cl))
    cl.mat <- purrr::map(cl.mat, cbind, cons.cl)
    an <- c(an, colnames(cons.cl))
  }

  # Internal indices
  ind.int <- purrr::map(cl.mat, ~ {
    x <- as.matrix(data)
    data.frame(
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
      mutate_all(funs(structure(., names = an)))
  })

  # Graph all plotting functions
  if (plot) {
    graph_all(E)
  }

  # Only calculate external indices if a reference is given
  if (!is.null(ref.cl)) {
    mat.ext <- cl.mat %>%
      magrittr::extract(match(k, names(.)))
    ind.ext <- purrr::map(mat.ext, ~ data.frame(
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

  # Only trim if specified and more than one algorithm
  if (dim(E)[3] > 1 & trim) {
    trim.obj <- purrr::map(k, ~ consensus_trim(E = E, ii = ind.int, k = .x,
                                               k.method = k.method,
                                               reweigh = reweigh, n = n)) %>%
      purrr::transpose() %>%
      purrr::map_at(c("alg.keep", "alg.remove"), ~ unlist(unique(.x)))
  } else {
    trim.obj <- list(alg.keep = an,
                     alg.remove = character(0),
                     rank.agg = list(NULL),
                     top.list = list(NULL),
                     data.new = list(E))
  }

  # Reorder ind.int (and ind.ext if not NULL) by top.list order if trimmed
  if (all(purrr::map_lgl(trim.obj$top.list, ~ !is.null(.x)))) {
    ind.int <- purrr::map2(ind.int, trim.obj$top.list,
                           ~ arrange(.x, match(.y, Algorithms)))
    if (!is.null(ind.ext))
      ind.ext <- purrr::map2(ind.ext, trim.obj$top.list,
                             ~ arrange(.x, match(.y, Algorithms)))
  }

  return(list(k = k, pac = pac, internal = ind.int, external = ind.ext,
              trim = trim.obj))
}

#' @param E consensus object from \code{consensus_evaluate}
#' @param ii internal indices object from \code{consensus_evaluate}
#' @param k chosen value(s) of k from \code{consensus_evaluate}
#' @noRd
consensus_trim <- function(E, ii, k, k.method, reweigh, n) {
  k <- as.character(k)
  zk <- ii[[k]]
  alg.all <- dimnames(E)[[3]]

  # Separate algorithms into those from clusterCrit (main), and (others)
  z.main <- zk %>%
    magrittr::extract(!names(.) %in% c("Algorithms", "Compactness",
                                       "Connectivity") &
              purrr::map_lgl(., ~ all(!is.nan(.x))))
  z.other <- zk %>%
    magrittr::extract(c("Compactness", "Connectivity"))

  # Which algorithm is the best for each index?
  bests <- purrr::map2_int(z.main, names(z.main), clusterCrit::bestCriterion)

  max.bests <- z.main %>%
    magrittr::extract(purrr::map_int(., which.max) == bests) %>%
    cbind(z.other) %>%
    multiply_by(-1)
  min.bests <- z.main %>%
    magrittr::extract(purrr::map_int(., which.min) == bests)

  # Determine trimmed ensemble using rank aggregation, only if there are more
  # algorithms than we want to keep
  if (length(alg.all) <= n) {
    rank.agg <- NULL
    top.list <- NULL
    alg.keep <- alg.all
  } else {
    rank.agg <- cbind(max.bests, min.bests) %>%
      scale(center = FALSE, scale = TRUE) %>%
      as.data.frame() %>%
      purrr::map_df(~ alg.all[order(.x, sample(.x))]) %>%
      t()
    top.list <- rank.agg %>%
      RankAggreg::RankAggreg(., ncol(.), method = "GA", verbose = FALSE) %>%
      use_series("top.list")
    alg.keep <- top.list[seq_len(n)]
  }
  alg.remove <- as.character(alg.all[!(alg.all %in% alg.keep)])
  E.trim <- E[, , alg.keep, k, drop = FALSE]

  # Reweigh only if specified, more than 1 algorithm is kept, trimming done
  if (reweigh && length(alg.keep) > 1 && !is.null(top.list)) {

    # Filter after knowing which to keep
    ak <- match(alg.keep, alg.all)
    max.bests <- max.bests[ak, ]
    min.bests <- z.main %>%
      magrittr::extract(ak, purrr::map_int(., which.min) == bests) %>%
      purrr::map_df(~ sum(.x) - .x)

    # Create multiples of each algorithm proportion to weight
    # Divide multiples by greatest common divisor to minimize number of copies
    multiples <- cbind(as.matrix(max.bests), as.matrix(min.bests)) %>%
      prop.table(2) %>%
      rowMeans() %>%
      multiply_by(100) %>%
      round(0) %>%
      divide_by(Reduce("gcd", .)) %>%
      set_names(alg.keep)

    # Generate multiples for each algorithm, adding back dimnames metadata
    E.trim <- purrr::array_branch(E.trim, c(3, 4)) %>%
      purrr::map2(., multiples, ~ rep(list(.x), .y)) %>%
      purrr::map(abind::abind, along = 3) %>%
      abind::abind(along = 3) %>%
      abind::abind(along = 4)
    dimnames(E.trim) <-
      list(NULL,
           dimnames(E.trim)[[2]],
           purrr::map2(names(multiples), multiples, rep) %>%
             purrr::flatten_chr(),
           k)
  }

  # If k.method is to select "all", need to add suffixes to algorithms
  if (!is.null(k.method) && k.method == "all") {
    alg.keep <- paste0(alg.keep, " k=", k)
    if (length(alg.remove) > 0) alg.remove <- paste0(alg.remove, " k=", k)
    dimnames(E.trim)[[3]] <- paste0(dimnames(E.trim)[[3]], " k=", k)
  }
  return(list(alg.keep = alg.keep,
              alg.remove = alg.remove,
              rank.agg = rank.agg,
              top.list = top.list,
              data.new = E.trim))
}

#' Recursively find the greater common divisor of two numbers
#' @noRd
gcd <- function(x, y) {
  r <- x %% y
  return(ifelse(r, gcd(y, r), y))
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

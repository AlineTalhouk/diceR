#' Evaluate, trim, and reweigh algorithms
#'
#' Evaluates algorithms on internal/external validation indices. Poor performing
#' algorithms can be trimmed from the ensemble. The remaining algorithms can be
#' given weights before use in consensus functions.
#'
#' This function always returns internal indices. If `ref.cl` is not `NULL`,
#' external indices are additionally shown. Relevant graphical displays are also
#' outputted. Algorithms are ranked across internal indices using Rank
#' Aggregation. Only the top `n` algorithms are kept, the rest are trimmed.
#'
#' @param data data matrix with rows as samples and columns as variables
#' @param cons.cl matrix of cluster assignments from consensus functions such as
#'   `kmodes` and `majority_voting`
#' @param ref.cl reference class
#' @param k.method determines the method to choose k when no reference class is
#'   given. When `ref.cl` is not `NULL`, k is the number of distinct classes of
#'   `ref.cl`. Otherwise the input from `k.method` chooses k. The default is to
#'   use the PAC to choose the best k(s). Specifying an integer as a
#'   user-desired k will override the best k chosen by PAC. Finally, specifying
#'   "all" will produce consensus results for all k. The "all" method is
#'   implicitly performed when there is only one k used.
#' @param plot logical; if `TRUE`, `graph_all` is called
#' @param trim logical; if `TRUE`, algorithms that score low on internal indices
#'   will be trimmed out
#' @param reweigh logical; if `TRUE`, after trimming out poor performing
#'   algorithms, each algorithm is reweighed depending on its internal indices.
#' @param n an integer specifying the top `n` algorithms to keep after trimming
#'   off the poor performing ones using Rank Aggregation. If the total number of
#'   algorithms is less than `n` no trimming is done.
#' @inheritParams consensus_combine
#' @return `consensus_evaluate` returns a list with the following elements
#' * `k`: if `ref.cl` is not `NULL`, this is the number of distinct classes
#'   in the reference; otherwise the chosen `k` is determined by the one giving
#'   the largest mean PAC across algorithms
#' * `pac`: a data frame showing the PAC for each combination of algorithm
#'   and cluster size
#' * `ii`: a list of data frames for all k showing internal evaluation
#'   indices
#' * `ei`: a data frame showing external evaluation indices for `k`
#' * `trim.obj`: A list with 4 elements
#'      * `alg.keep`: algorithms kept
#'      * `alg.remove`: algorithms removed
#'      * `rank.matrix`: a matrix of ranked algorithms for every internal
#'        evaluation index
#'      * `top.list`: final order of ranked algorithms
#'      * `E.new`: A new version of a `consensus_cluster` data object
#' @export
#' @examples
#' # Consensus clustering for multiple algorithms
#' set.seed(911)
#' x <- matrix(rnorm(500), ncol = 10)
#' CC <- consensus_cluster(x, nk = 3:4, reps = 10, algorithms = c("ap", "km"),
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
  cl.mat <- consensus_combine(E, element = "class") %>%
    purrr::map(as.data.frame)
  cons.mat <- consensus_combine(E, element = "matrix")

  # Calculate PAC and choose k
  pac <- cons.mat %>%
    purrr::modify_depth(2, PAC) %>%
    purrr::map_df(data.frame, .id = "k")
  k <- choose_k(ref.cl, k.method, pac)

  # If matrix of cluster assignments from cons.funs given, cbind to cl.mat
  an <- dimnames(E)[[3]]
  if (!is.null(cons.cl)) {
    assertthat::assert_that(is.matrix(cons.cl))
    cl.mat <- purrr::map(cl.mat, cbind, cons.cl)
    an <- c(an, colnames(cons.cl))
  }

  # Internal indices
  ii <- cl.mat %>% purrr::map(ivi_table, data = data)

  # External indices if reference is given
  ei <- ref.cl %>%
    purrr::when(
      !is.null(.) ~ cl.mat %>%
        magrittr::extract(match(k, names(.))) %>% # for chosen k only
        purrr::map(evi_table, ref.cl = ref.cl),
      TRUE ~ NULL
    )

  # Only trim if more than one algorithm and trim is specified
  if (dim(E)[3] > 1 & trim) {
    trim.obj <- purrr::map(k, consensus_trim, E = E, ii = ii,
                           k.method = k.method, reweigh = reweigh, n = n) %>%
      purrr::transpose() %>%
      purrr::map_at(c("alg.keep", "alg.remove"),
                    ~ unlist(unique(.x)))
  } else {
    trim.obj <- list(alg.keep = an,
                     alg.remove = character(0),
                     rank.matrix = list(NULL),
                     top.list = list(NULL),
                     E.new = list(E))
  }

  # Reorder ii (and ei if not NULL) by top.list order if trimmed
  if (all(purrr::map_lgl(trim.obj$top.list, ~ !is.null(.x)))) {
    ii <- purrr::map2(ii, trim.obj$top.list,
                      ~ dplyr::arrange(.x, match(.y, Algorithms)))
    if (!is.null(ei))
      ei <- purrr::map2(ei, trim.obj$top.list,
                        ~ dplyr::arrange(.x, match(.y, Algorithms)))
  }

  if (plot) graph_all(E) # Graph all plotting functions

  dplyr::lst(k, pac, ii, ei, trim.obj)
}

#' @param E consensus object from `consensus_evaluate`
#' @param ii internal indices object from `consensus_evaluate`
#' @param k chosen value(s) of k from `consensus_evaluate`
#' @noRd
consensus_trim <- function(E, ii, k, k.method, reweigh, n) {
  # Extract ii only for chosen k
  k <- as.character(k)
  ii <- ii[[k]]
  alg.all <- as.character(ii$Algorithms)

  # Rank algorithms on internal indices, store rank matrix and top alg list
  rank.obj <- consensus_rank(ii, n)
  rank.matrix <- rank.obj$rank.matrix
  top.list <- rank.obj$top.list

  # Algorithms to keep and remove based on top-ranked list
  alg.keep <- top.list %>% purrr::when(is.null(.) ~ alg.all,
                                       TRUE ~ .[seq_len(n)])
  alg.remove <- as.character(alg.all[!(alg.all %in% alg.keep)])
  E.new <- E[, , alg.keep, k, drop = FALSE]

  # Reweigh only if specified, more than 1 algorithm is kept, trimming done
  if (reweigh && length(alg.keep) > 1 && !is.null(top.list)) {
    E.new <- consensus_reweigh(E.new, rank.obj, alg.keep, alg.all)
  }

  # If k.method is to select "all", need to add suffixes to algorithms
  if (!is.null(k.method) && k.method == "all") {
    alg.keep <- paste_k(alg.keep, k)
    alg.remove <- alg.remove %>%
      purrr::when(length(.) > 0 ~ paste_k(., k), TRUE ~ .)
    dimnames(E.new)[[3]] <- paste_k(dimnames(E.new)[[3]], k)
  }
  dplyr::lst(alg.keep, alg.remove, rank.matrix, top.list, E.new)
}

#' Rank based on internal validity indices
#' @noRd
consensus_rank <- function(ii, n) {
  # Extract internal indices from clusterCrit and remove NaN idx
  ii.cc <- ii %>%
    magrittr::extract(
      !names(.) %in% c("Algorithms", "Compactness", "Connectivity") &
        purrr::map_lgl(., ~ all(!is.nan(.x)))
    )

  # Which algorithm is the best for each index?
  bests <- purrr::imap_int(ii.cc, clusterCrit::bestCriterion)
  max.bests <- ii.cc %>%
    magrittr::extract(purrr::map_int(., which.max) == bests) %>%
    magrittr::multiply_by(-1)
  min.bests <- ii.cc %>%
    magrittr::extract(purrr::map_int(., which.min) == bests) %>%
    cbind(ii[c("Compactness", "Connectivity")]) # these two need to be minimized

  # Determine trimmed ensemble using rank aggregation
  if (nrow(ii) <= n) {
    rank.matrix <- top.list <- NULL
  } else {
    rank.matrix <- cbind(max.bests, min.bests) %>%
      scale(center = FALSE, scale = TRUE) %>%
      as.data.frame() %>%
      purrr::map_dfc(~ ii$Algorithms[order(.x, sample(length(.x)))]) %>%
      t()
    top.list <- RankAggreg::RankAggreg(rank.matrix, ncol(rank.matrix),
                                       method = "GA", verbose = FALSE)$top.list
  }
  dplyr::lst(max.bests, min.bests, rank.matrix, top.list)
}

#' Reweigh the algorithms in the ensemble if some were trimmed out
#' @noRd
consensus_reweigh <- function(E.new, rank.obj, alg.keep, alg.all) {
  # Create multiples of each algorithm proportion to weight
  multiples <- rank.obj %>%
    magrittr::extract(c("max.bests", "min.bests")) %>%
    dplyr::bind_cols() %>% # Recombine internal validity indices
    dplyr::select_if(~ !any(is.infinite(.))) %>% # Remove indices with Inf
    magrittr::extract(match(alg.keep, alg.all), ) %>% # Extract algs to keep
    as.matrix() %>% # Cannot calculate proportions on data.frame
    prop.table(2) %>% # Proportion for each index
    rowMeans() %>% # Average proportion across indices
    magrittr::multiply_by(100) %>% # Percentage
    round(0) %>% # Round to closest whole number
    magrittr::divide_by(Reduce(`gcd`, .)) %>% # Divide by gcd to minimize copies
    purrr::set_names(alg.keep) # reassign names of algs to keep

  # Generate multiples for each algorithm, updating dimnames 3rd dimension
  E.new %>%
    purrr::array_branch(c(3, 4)) %>%
    purrr::map2(multiples, ~ rep(list(.x), .y)) %>%
    purrr::map(abind::abind, along = 3) %>%
    abind::abind(along = 3) %>%
    abind::abind(along = 4) %>%
    `dimnames<-`(append(dimnames(E.new)[-3],
                        list(purrr::flatten_chr(
                          purrr::imap(multiples, ~ rep(.y, .x))
                        )),
                        2))
}

#' Table of internal validity indices for each algorithm
#' @param cl.df data frame of cluster assignments for each algorithm
#' @param data data frame with rows as samples and columns as variables
#' @noRd
ivi_table <- function(cl.df, data) {
  ndata <- apply(data, 2, function(x) as.numeric(as.character(x)))
  data.frame(
    Algorithms = colnames(cl.df),
    cl.df %>% purrr::map_df(
      clusterCrit::intCriteria,
      traj = ndata,
      crit = c("Calinski_Harabasz", "Dunn", "PBM", "Tau", "Gamma", "C_index",
               "Davies_Bouldin", "McClain_Rao", "SD_Dis", "Ray_Turi", "G_plus",
               "Silhouette", "S_Dbw")),
    Compactness = cl.df %>% purrr::map_dbl(compactness, data = data),
    Connectivity = cl.df %>% purrr::map_dbl(
      ~ clValid::connectivity(Data = ndata, clusters = .))
  ) %>%
    dplyr::mutate_all(list(~ structure(., names = colnames(cl.df))))
}

#' Table of external validity indices for each algorithm
#' @param cl.df data frame of cluster assignments for each algorithm
#' @param ref.cl reference class
#' @noRd
evi_table <- function(cl.df, ref.cl) {
  data.frame(
    Algorithms = colnames(cl.df),
    cl.df %>% purrr::map_df(
      clusterCrit::extCriteria,
      part2 = ref.cl,
      crit = c("Hubert", "Jaccard", "McNemar", "Rand")
    ),
    NMI = cl.df %>% purrr::map_dbl(ev_nmi, ref.lab = ref.cl)
  ) %>%
    dplyr::inner_join(
      cl.df %>%
        purrr::map_dfr(ev_confmat, ref.lab = ref.cl, .id = "Algorithms") %>%
        dplyr::select(-.data$.estimator) %>%
        tidyr::pivot_wider(names_from = ".metric", values_from = ".estimate"),
      by = "Algorithms"
    ) %>%
    dplyr::mutate_all(list(~ structure(., names = colnames(cl.df))))
}

#' Choose k using PAC
#' @noRd
choose_k <- function(ref.cl, k.method, pac) {
  # If reference given, k is number of distinct classes
  if (!is.null(ref.cl)) {
    k <- dplyr::n_distinct(ref.cl)
    # Otherwise k is chosen using the following methods
  } else if (is.null(k.method)) {
    k <- pac %>%
      magrittr::use_series("k") %>%
      magrittr::extract(pac[, -1, drop = FALSE] %>%
                          apply(2, which.min) %>%
                          unlist() %>%
                          table() %>%
                          magrittr::extract(magrittr::is_in(., max(.))) %>%
                          names() %>%
                          as.numeric()) %>%
      as.integer() %>%
      min()
  } else if (any(nrow(pac) == 1 | k.method == "all")) {
    k <- pac$k
  } else if (purrr::is_scalar_integer(k.method) ||
             purrr::is_scalar_double(k.method)) {
    k <- k.method
  } else {
    stop("Invalid input. Check documentation for possible options.")
  }
  as.integer(k)
}

#' Add cluster size suffix to algorithms when k.method == "all"
#' @noRd
paste_k <- function(a, k) {
  paste0(a, " k=", k)
}

#' Recursively find the greater common divisor of two numbers
#' @noRd
gcd <- function(x, y) {
  r <- x %% y
  ifelse(r, gcd(y, r), y)
}

#' Compactness Measure
#'
#' Compute the compactness validity index for a clustering result.
#'
#' This index is agnostic to any reference clustering results, calculating
#' cluster performance on the basis of compactness and separability. Smaller
#' values indicate a better clustering structure.
#'
#' @param data a dataset with rows as observations, columns as variables
#' @param labels a vector of cluster labels from a clustering result
#' @return the compactness score
#' @author Derek Chiu
#' @references MATLAB function `valid_compactness` by Simon Garrett in
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
  cp / length(labels)
}

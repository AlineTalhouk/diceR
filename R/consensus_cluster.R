#' Consensus clustering
#'
#' Runs consensus clustering across subsamples of the data, clustering
#' algorithms, and cluster sizes.
#'
#' See examples for how to use custom algorithms and distance functions. The
#' default clustering algorithms provided are:
#' \itemize{
#'   \item{"nmf": }{Nonnegative Matrix Factorization (using Kullback-Leibler
#'   Divergence or Euclidean distance; See Note for specifications.)}
#'   \item{"hc": }{Hierarchical Clustering}
#'   \item{"diana": }{DIvisive ANAlysis Clustering}
#'   \item{"km": }{K-Means Clustering}
#'   \item{"pam": }{Partition Around Medoids}
#'   \item{"ap": }{Affinity Propagation}
#'   \item{"sc": }{Spectral Clustering using Radial-Basis kernel function}
#'   \item{"gmm": }{Gaussian Mixture Model using Bayesian Information Criterion
#'   on EM algorithm}
#'   \item{"block": }{Biclustering using a latent block model}
#'   \item{"som": }{Self-Organizing Map (SOM) with Hierarchical Clustering}
#'   \item{"cmeans": }{Fuzzy C-Means Clustering}
#'   \item{"hdbscan": }{Hierarchical Density-based Spatial Clustering of
#'   Applications with Noise (HDBSCAN)}
#' }
#'
#' The progress bar increments on every unit of \code{reps}.
#'
#' @note The \code{nmf.method} defaults are "brunet" (Kullback-Leibler
#'   divergence) and "lee" (Euclidean distance). When "hdbscan" is chosen as an
#'   algorithm to use, its results are excluded from the rest of the consensus
#'   clusters. This is because there is no guarantee that the cluster assignment
#'   will have every sample clustered; more often than not there will be noise
#'   points or outliers. In addition, the number of distinct clusters may not
#'   even be equal to \code{nk}.
#'
#' @param data data matrix with rows as samples and columns as variables
#' @param nk number of clusters (k) requested; can specify a single integer or a
#'   range of integers to compute multiple k
#' @param p.item proportion of items to be used in subsampling within an
#'   algorithm
#' @param reps number of subsamples
#' @param algorithms vector of clustering algorithms for performing consensus
#'   clustering. Must be any number of the following: "nmf", "hc", "diana",
#'   "km", "pam", "ap", "sc", "gmm", "block", "som", "cmeans", "hdbscan". A
#'   custom clustering algorithm can be used.
#' @param nmf.method specify NMF-based algorithms to run. By default the
#'   "brunet" and "lee" algorithms are called. See \code{\link[NMF]{nmf}} for
#'   details.
#' @param xdim x dimension of the SOM grid
#' @param ydim y dimension of the SOM grid
#' @param rlen the number of times the complete data set will be presented to
#'   the SOM network.
#' @param alpha SOM learning rate, a vector of two numbers indicating the amount
#'   of change. Default is to decline linearly from 0.05 to 0.01 over
#'   \code{rlen} updates. Not used for the batch algorithm.
#' @param minPts minimum size of clusters for HDBSCAN. Default is 5.
#' @param distance a vector of distance functions. Defaults to "euclidean".
#'   Other options are given in \code{\link[stats]{dist}}. A custom distance
#'   function can be used.
#' @param prep.data Prepare the data on the "full" dataset, the
#'   "sampled" dataset, or "none" (default).
#' @inheritParams prepare_data
#' @param progress logical; should a progress bar be displayed?
#' @param seed.nmf random seed to use for NMF-based algorithms
#' @param seed.data seed to use to ensure each algorithm operates on the same
#'   set of subsamples
#' @param file.name if not \code{NULL}, the returned array will be saved at each
#'   iteration as well as at the end of the function call to an \code{rds}
#'   object with \code{file.name} as the file name.
#' @param time.saved logical; if \code{TRUE}, the date saved is appended to
#'   \code{file.name}. Only applicable when \code{file.name} is not \code{NULL}.
#' @return An array of dimension \code{nrow(x)} by \code{reps} by
#'   \code{length(algorithms)} by \code{length(nk)}. Each cube of the array
#'   represents a different k. Each slice of a cube is a matrix showing
#'   consensus clustering results for algorithms. The matrices have a row for
#'   each sample, and a column for each subsample. Each entry represents a class
#'   membership.
#'
#'   When "hdbscan" is part of \code{algorithms}, we do not include its clustering
#'   array in the consensus result. Instead, we report two summary statistics
#'   as attributes: the proportion of outliers and the number of clusters.
#' @author Derek Chiu, Aline Talhouk
#' @importFrom mclust mclustBIC
#' @export
#' @examples
#' data(hgsc)
#' dat <- hgsc[1:100, 1:50]
#'
#' # Custom distance function
#' manh <- function(x) {
#'   stats::dist(x, method = "manhattan")
#' }
#'
#' # Custom clustering algorithm
#' agnes <- function(d, k) {
#'   return(as.integer(stats::cutree(cluster::agnes(d, diss = TRUE), k)))
#' }
#'
#' assign("agnes", agnes, 1)
#'
#' cc <- consensus_cluster(dat, reps = 6, algorithms = c("pam", "agnes"),
#' distance = c("euclidean", "manh"), progress = FALSE)
#' str(cc)
consensus_cluster <- function(data, nk = 2:4, p.item = 0.8, reps = 1000,
                              algorithms = NULL,
                              nmf.method = c("brunet", "lee"),
                              xdim = NULL, ydim = NULL, rlen = 200,
                              alpha = c(0.05, 0.01), minPts = 5,
                              distance = "euclidean",
                              prep.data = c("none", "full", "sampled"),
                              scale = TRUE, type = c("conventional", "robust"),
                              min.var = 1, progress = TRUE,
                              seed.nmf = 123456, seed.data = 1,
                              file.name = NULL, time.saved = FALSE) {
  prep.data <- match.arg(prep.data)
  if (prep.data == "full")
    data <- prepare_data(data, scale = scale, type = type, min.var = min.var)
  algorithms <- algorithms %||% ALG_NAMES  # Use all if none are specified

  # Store consensus dimensions for calculating progress bar increments/offsets
  lnk <- length(nk)
  algs <- dplyr::lst(NALG, DALG, OALG) %>%
    purrr::map(~ algorithms[algorithms %in% .x])
  lnmf <- length(algs$NALG) * length(nmf.method)
  ldist <- length(algs$DALG) * length(distance)
  lother <- length(algs$OALG)
  if (progress) {
    pb <- utils::txtProgressBar(max = lnk * (lnmf + ldist + lother) * reps,
                                style = 3)
  } else {
    pb <- NULL
  }

  # Argument lists: Common, NMF, Distance, Other
  cargs <- dplyr::lst(data, nk, p.item, reps, seed.data, prep.data, scale, type,
                      min.var, progress, pb)
  nargs <- c(cargs, dplyr::lst(algs = algs$NALG, nmf.method, seed.nmf))
  dargs <- c(cargs, dplyr::lst(algs = algs$DALG, distance,
                               offset = lnk * lnmf * reps))
  oargs <- c(cargs, dplyr::lst(algs = algs$OALG, xdim, ydim, rlen, alpha,
                               minPts, offset = lnk * (lnmf + ldist) * reps))
  args <- list(nargs, dargs, oargs)

  # Run cc on all algorithms, combine on 3rd dim, HDBSCAN manipulation
  fun <- list(cc_nmf, cc_dist, cc_other)
  arr_all <- purrr::pmap(list(fun, args), cc) %>%
    abind::abind(along = 3) %>%
    hdbscan_summarize(algorithms)

  if (!is.null(file.name)) {
    if (time.saved) {
      file.name <- paste0(file.name, "_",
                          format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
    }
    saveRDS(arr_all, file = paste0(file.name, ".rds"))
  }
  arr_all
}

#' If algs from each group exist, invoke respective cc fun, otherwise NULL
#' @noRd
cc <- function(fun, args) {
  length(args$algs) %>% purrr::when(. > 0 ~ purrr::invoke(fun, args), ~ NULL)
}

#' Cluster NMF-based algorithms
#' @noRd
cc_nmf <- function(data, nk, p.item, reps, algs, nmf.method, seed.nmf,
                   seed.data, prep.data, scale, type, min.var, progress,
                   pb) {
  x_nmf <- nmf_transform(data)
  n <- nrow(data)
  alg <- paste(toupper(algs), Hmisc::capitalize(nmf.method), sep = "_")
  arr_nmf <- init_array(data, reps, alg, nk)

  for (k in seq_along(nk)) {
    for (j in seq_along(nmf.method)) {
      set.seed(seed.data)
      for (i in seq_len(reps)) {
        ind.new <- sample(n, floor(n * p.item))
        # In case the subsample has all-zero vars, remove them to speed up comp
        x <- x_nmf[ind.new, ] %>% magrittr::extract(colSums(.) != 0)
        if (prep.data == "sampled") {
          x <- x %>%
            prepare_data(scale = scale, type = type, min.var = min.var) %>%
            nmf_transform()
        }
        arr_nmf[ind.new, i, j, k] <- nmf(x, nk[k], nmf.method[j], seed.nmf)

        if (progress) {
          value <- (k - 1) * length(nmf.method) * reps + (j - 1) * reps + i
          utils::setTxtProgressBar(pb, value)
        }
      }
    }
  }
  arr_nmf
}

#' Cluster algorithms with dissimilarity specification
#' @noRd
cc_dist <- function(data, nk, p.item, reps, algs, distance, seed.data,
                    prep.data, scale, type, min.var, progress, pb,
                    offset) {
  n <- nrow(data)
  alg <- apply(expand.grid(Hmisc::capitalize(distance),
                           toupper(algs)),
               1, function(x) paste0(x[2], "_", x[1]))
  arr_dist <- init_array(data, reps, alg, nk)

  for (k in seq_along(nk)) {
    for (j in seq_along(algs)) {
      for (d in seq_along(distance)) {
        set.seed(seed.data)
        for (i in seq_len(reps)) {
          ind.new <- sample(n, floor(n * p.item))
          x <- data[ind.new, ]

          if (prep.data == "sampled") {
            x <- prepare_data(x, scale = scale, type = type, min.var = min.var)
          }
          dists <- distances(x, distance[d])
          arr_dist[ind.new, i, (j - 1) * length(distance) + d, k] <-
            get(algs[j])(dists[[1]], nk[k])  # get() for custom funs
          if (progress) {
            value <- (k - 1) * length(algs) * length(distance) * reps +
              (j - 1) * length(distance) * reps + (d - 1) * reps + i + offset
            utils::setTxtProgressBar(pb, value)
          }
        }
      }
    }
  }
  arr_dist
}

#' Cluster other algorithms
#' @noRd
cc_other <- function(data, nk, p.item, reps, algs, xdim, ydim, rlen,
                     alpha, seed.data, prep.data, scale, type, min.var,
                     progress, pb, minPts, offset) {
  n <- nrow(data)
  alg <- toupper(algs)
  arr_other <- init_array(data, reps, alg, nk)

  for (k in seq_along(nk)) {
    for (j in seq_along(algs)) {
      set.seed(seed.data)
      for (i in seq_len(reps)) {
        ind.new <- sample(n, floor(n * p.item))
        x <- data[ind.new, ]

        if (prep.data == "sampled") {
          x <- prepare_data(x, scale = scale, type = type, min.var = min.var)
        }
        arr_other[ind.new, i, j, k] <-
          switch(algs[j],
                 ap = ap(x, nk[k]),
                 sc = sc(x, nk[k]),
                 gmm = gmm(x, nk[k]),
                 block = block(x, nk[k]),
                 som = som(x, nk[k], xdim, ydim, rlen, alpha),
                 cmeans = cmeans(x, nk[k]),
                 hdbscan = hdbscan(x, minPts)
          )
        if (progress) {
          value <- (k - 1) * length(algs) * reps + (j - 1) * reps + i + offset
          utils::setTxtProgressBar(pb, value)
        }
      }
    }
  }
  arr_other
}

#' Initialize array to store consensus clustering results
#' @noRd
init_array <- function(data, r, a, k) {
  rn <- rownames(data) %||% seq_len(nrow(data))
  dn <- list(rn, paste0("R", seq_len(r)), a, k)
  array(NA_integer_, dim = purrr::map_int(dn, length), dimnames = dn)
}

#' Return a list of distance matrices
#' @param x data matrix
#' @param dist a character vector of distance methods taken from stats::dist, or
#'   "spearman", or a custom distance function in the current environment
#' @noRd
distances <- function(x, dist) {
  # Change partial matching distance methods from stats::dist to full names
  M <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
  dist <- ifelse(dist %pin% M, M[pmatch(dist, M)], dist)

  # Check if spearman distance is used (starts with string)
  sp.idx <- purrr::map_lgl(paste0("^", dist), grepl, "spearman")
  if (any(sp.idx)) {
    spear <- purrr::set_names(list(spearman_dist(x)), "spearman")
    d <- dist[!sp.idx]
  } else {
    spear <- NULL
    d <- dist
  }

  # If only spearman requested
  if (length(d) == 0) {
    return(spear)
  } else {
    # Identify custom distance functions from those found in stats::dist
    check <- d %>%
      purrr::map(~ try(stats::dist(x = x, method = .x), silent = TRUE)) %>%
      purrr::set_names(d)
    is.error <- purrr::map_lgl(check, inherits, "try-error")
    succeeded <- which(!is.error)
    failed <- which(is.error)

    # Search for custom function in parent environments
    if (length(failed)) {
      custom <- check[failed] %>%
        names() %>%
        purrr::map(~ get(.x)(x)) %>%
        purrr::set_names(names(failed))
    } else {
      custom <- NULL
    }

    # Combine distances into list and return in same order as dist argument
    dlist <- c(spear, check[succeeded], custom) %>%
      magrittr::extract(pmatch(dist, names(.)))
    return(dlist)
  }
}

#' Calculate pairwise Spearman correlational distances using
#' bioDist::spearman.dist defaults
#' @references
#' https://github.com/Bioconductor-mirror/bioDist/blob/master/R/spearman.dist.R
#'
#' @noRd
spearman_dist <- function(x) {
  x %>%
    t() %>%
    stats::cor(method = "spearman") %>%
    abs() %>%
    magrittr::subtract(1, .) %>%
    magrittr::extract(lower.tri(.)) %>%
    `attributes<-`(
      list(
        Size = nrow(x),
        Labels = rownames(x),
        Diag = FALSE,
        Upper = FALSE,
        methods = "spearman",
        class = "dist"
      )
    )
}

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
                              scale = TRUE,
                              type = c("conventional", "robust", "tsne",
                                       "largevis"),
                              min.var = 1, progress = TRUE,
                              seed.nmf = 123456, seed.data = 1,
                              file.name = NULL, time.saved = FALSE) {
  prep.data <- match.arg(prep.data)
  if (prep.data == "full")
    data <- prepare_data(data, scale = scale, type = type, min.var = min.var)
  algorithms <- algorithms %||% ALG_NAMES  # Use all if none are specified

  # Calculate total number of algorithms, including custom ones into DALG
  calgs <- algorithms[!algorithms %in% ALG_NAMES]
  algs <- dplyr::lst(NALG, DALG = c(DALG, calgs), OALG) %>%
    purrr::map(~ algorithms[algorithms %in% .x])
  lalg <- lengths(algs) * lengths(list(nmf.method, distance, 1))
  n <- nrow(data)

  if (progress) {
    pb <- progress::progress_bar$new(
      format = "Clustering Algorithm :num of :den: :alg (k = :k) [:bar] :percent eta: :eta",
      total = length(nk) * sum(lalg) * reps,
      clear = FALSE
    )
  } else {
    pb <- NULL
  }

  # Argument lists: Common, NMF, Distance, Other
  cargs <- dplyr::lst(data, nk, p.item, reps, seed.data, prep.data, scale, type,
                      min.var, pb, lalg, n)
  nargs <- dplyr::lst(algs = algs$NALG, nmf.method, seed.nmf)
  dargs <- dplyr::lst(algs = algs$DALG, distance)
  oargs <- dplyr::lst(algs = algs$OALG, xdim, ydim, rlen, alpha, minPts)
  args <- purrr::map(list(nargs, dargs, oargs), ~ c(cargs, .))

  # Run cc on all algorithms, combine on 3rd dim
  fun <- list(cc_nmf, cc_dist, cc_other)
  arr_all <- purrr::pmap(list(fun, args), cc) %>% abind::abind(along = 3)
  if ("hdbscan" %in% algorithms) {
    arr_all <- hdbscan_summarize(arr_all)  # HDBSCAN summaries
  }

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
                   seed.data, prep.data, scale, type, min.var, pb, lalg, n) {
  alg <- paste(toupper(algs), Hmisc::capitalize(nmf.method), sep = "_")
  arr_nmf <- init_array(data, reps, alg, nk)
  x_nmf <- nmf_transform(data)

  for (j in seq_along(nmf.method)) {
    for (k in seq_along(nk)) {
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
        if (!is.null(pb)) {
          pb$tick(tokens = list(num = j, den = sum(lalg), alg = alg[j],
                                k = nk[k]))
        }
        arr_nmf[ind.new, i, j, k] <- nmf(x, nk[k], nmf.method[j], seed.nmf)
      }
    }
  }
  arr_nmf
}

#' Cluster algorithms with dissimilarity specification
#' @noRd
cc_dist <- function(data, nk, p.item, reps, algs, distance, seed.data,
                    prep.data, scale, type, min.var, pb, lalg, n) {
  alg <- paste(rep(toupper(algs), each = length(distance)),
               rep(Hmisc::capitalize(distance), length(algs)),
               sep = "_")
  arr_dist <- init_array(data, reps, alg, nk)

  for (j in seq_along(algs)) {
    for (k in seq_along(nk)) {
      for (d in seq_along(distance)) {
        set.seed(seed.data)
        for (i in seq_len(reps)) {
          ind.new <- sample(n, floor(n * p.item))
          x <- data[ind.new, ]
          if (prep.data == "sampled") {
            x <- prepare_data(x, scale = scale, type = type, min.var = min.var)
          }
          dists <- cdist(x, distance[d])
          a <- (j - 1) * length(distance) + d
          if (!is.null(pb)) {
            pb$tick(tokens = list(num = j + lalg["NALG"], den = sum(lalg),
                                  alg = alg[a], k = nk[k]))
          }
          arr_dist[ind.new, i, a, k] <- get(algs[j])(dists, nk[k]) # custom
        }
      }
    }
  }
  arr_dist
}

#' Cluster other algorithms
#' @noRd
cc_other <- function(data, nk, p.item, reps, algs, xdim, ydim, rlen, alpha,
                     minPts, seed.data, prep.data, scale, type, min.var, pb,
                     lalg, n) {
  alg <- toupper(algs)
  arr_other <- init_array(data, reps, alg, nk)

  for (j in seq_along(algs)) {
    for (k in seq_along(nk)) {
      set.seed(seed.data)
      for (i in seq_len(reps)) {
        ind.new <- sample(n, floor(n * p.item))
        x <- data[ind.new, ]
        if (prep.data == "sampled") {
          x <- prepare_data(x, scale = scale, type = type, min.var = min.var)
        }
        if (!is.null(pb)) {
          pb$tick(tokens = list(num = j + sum(lalg[c("NALG", "DALG")]),
                                den = sum(lalg), alg = alg[j], k = nk[k]))
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

#' Return a clustering distance matrix object
#' @param x data matrix
#' @param dist a character string of distance methods taken from stats::dist, or
#'   "spearman", or a custom distance function in the current environment
#' @noRd
cdist <- function(x, dist) {
  # Change partial matching distance methods from stats::dist to full names
  M <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
  dist <- ifelse(dist %pin% M, M[pmatch(dist, M)], dist)

  # If spearman is specified, use spearman_dist()
  if (dist %pin% "spearman") return(spearman_dist(x))

  # Identify custom distance functions from those found in stats::dist
  d <- try(stats::dist(x = x, method = dist), silent = TRUE)

  # Run custom function in parent environments if not found
  if (inherits(d, "try-error")) get(dist)(x) else d
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

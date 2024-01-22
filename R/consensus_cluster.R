#' Consensus clustering
#'
#' Runs consensus clustering across subsamples of the data, clustering
#' algorithms, and cluster sizes.
#'
#' See examples for how to use custom algorithms and distance functions. The
#' default clustering algorithms provided are:
#' * "nmf": Nonnegative Matrix Factorization (using Kullback-Leibler Divergence
#' or Euclidean distance; See Note for specifications.)
#' * "hc": Hierarchical Clustering
#' * "diana": DIvisive ANAlysis Clustering
#' * "km": K-Means Clustering
#' * "pam": Partition Around Medoids
#' * "ap": Affinity Propagation
#' * "sc": Spectral Clustering using Radial-Basis kernel function
#' * "gmm": Gaussian Mixture Model using Bayesian Information Criterion on EM
#'   algorithm
#' * "block": Biclustering using a latent block model
#' * "som": Self-Organizing Map (SOM) with Hierarchical Clustering
#' * "cmeans": Fuzzy C-Means Clustering
#' * "hdbscan": Hierarchical Density-based Spatial Clustering of Applications
#'   with Noise (HDBSCAN)
#'
#' The progress bar increments on every unit of `reps`.
#'
#' @note The `nmf.method` options are "brunet" (Kullback-Leibler Divergence) and
#'   "lee" (Euclidean distance). When "hdbscan" is chosen as an algorithm to
#'   use, its results are excluded from the rest of the consensus clusters. This
#'   is because there is no guarantee that the cluster assignment will have
#'   every sample clustered; more often than not there will be noise points or
#'   outliers. In addition, the number of distinct clusters may not even be
#'   equal to `nk`.
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
#'   "brunet" and "lee" algorithms are called. See [NMF::nmf()] for details.
#' @param hc.method agglomeration method for hierarchical clustering. The
#'   the "average" method is used by default. See[stats::hclust()] for details.
#' @param xdim x dimension of the SOM grid
#' @param ydim y dimension of the SOM grid
#' @param rlen the number of times the complete data set will be presented to
#'   the SOM network.
#' @param alpha SOM learning rate, a vector of two numbers indicating the amount
#'   of change. Default is to decline linearly from 0.05 to 0.01 over `rlen`
#'   updates. Not used for the batch algorithm.
#' @param minPts minimum size of clusters for HDBSCAN. Default is 5.
#' @param distance a vector of distance functions. Defaults to "euclidean".
#'   Other options are given in [stats::dist()]. A custom distance function can
#'   be used.
#' @param abs only used for `distance = c("spearman", "pearson")`. If `TRUE`,
#'   the absolute value is first applied to the distance before subtracting from
#'   1, e.g., we use 1 - |SCD| instead of 1 - SCD for the spearman correlation
#'   distance.
#' @param prep.data Prepare the data on the "full" dataset, the "sampled"
#'   dataset, or "none" (default).
#' @inheritParams prepare_data
#' @param progress logical; should a progress bar be displayed?
#' @param seed.nmf random seed to use for NMF-based algorithms
#' @param seed.data seed to use to ensure each algorithm operates on the same
#'   set of subsamples
#' @param file.name if not `NULL`, the returned array will be saved at each
#'   iteration as well as at the end of the function call to an `rds` object
#'   with `file.name` as the file name.
#' @param time.saved logical; if `TRUE`, the date saved is appended to
#'   `file.name`. Only applicable when `file.name` is not `NULL`.
#' @return An array of dimension `nrow(x)` by `reps` by `length(algorithms)` by
#'   `length(nk)`. Each cube of the array represents a different k. Each slice
#'   of a cube is a matrix showing consensus clustering results for algorithms.
#'   The matrices have a row for each sample, and a column for each subsample.
#'   Each entry represents a class membership.
#'
#'   When "hdbscan" is part of `algorithms`, we do not include its clustering
#'   array in the consensus result. Instead, we report two summary statistics as
#'   attributes: the proportion of outliers and the number of clusters.
#' @author Derek Chiu, Aline Talhouk
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
                              hc.method = "average",
                              xdim = NULL, ydim = NULL, rlen = 200,
                              alpha = c(0.05, 0.01), minPts = 5,
                              distance = "euclidean",
                              abs = TRUE,
                              prep.data = c("none", "full", "sampled"),
                              scale = TRUE,
                              type = c("conventional", "robust", "tsne"),
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
    if (!requireNamespace("progress", quietly = TRUE)) {
      stop("Package \"progress\" is needed. Please install it.",
           call. = FALSE)
    } else {
      pb <- progress::progress_bar$new(
        format = "Clustering Algorithm :num of :den: :alg (k = :k) [:bar] :percent eta: :eta",
        total = length(nk) * sum(lalg) * reps,
        clear = FALSE
      )
    }
  } else {
    pb <- NULL
  }

  # Argument lists: Common, NMF, Distance, Other
  cargs <- dplyr::lst(data, nk, p.item, reps, seed.data, prep.data, scale, type,
                      min.var, pb, lalg, n)
  nargs <- dplyr::lst(algs = algs$NALG, nmf.method, seed.nmf)
  dargs <- dplyr::lst(algs = algs$DALG, distance, abs, hc.method)
  oargs <- dplyr::lst(algs = algs$OALG, xdim, ydim, rlen, alpha, minPts,
                      hc.method)
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
    saveRDS(arr_all, file = paste0(file.name, ".rds"), version = 2)
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
  alg <- paste(toupper(algs), stringr::str_to_title(nmf.method), sep = "_")
  arr <- init_array(data, reps, alg, nk)
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
        arr[ind.new, i, j, k] <- nmf(x, nk[k], nmf.method[j], seed.nmf)
      }
    }
  }
  arr
}

#' Cluster algorithms with dissimilarity specification
#' @noRd
cc_dist <- function(data, nk, p.item, reps, algs, distance, abs, hc.method,
                    seed.data, prep.data, scale, type, min.var, pb, lalg, n) {
  alg <- paste(rep(toupper(algs), each = length(distance)),
               rep(stringr::str_to_title(distance), length(algs)),
               sep = "_")
  arr <- init_array(data, reps, alg, nk)

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
          dists <- cdist(x = x, dist = distance[d], abs = abs)
          a <- (j - 1) * length(distance) + d
          if (!is.null(pb)) {
            pb$tick(tokens = list(num = j + lalg["NALG"], den = sum(lalg),
                                  alg = alg[a], k = nk[k]))
          }
          if (algs[j] == "hc") {
            arr[ind.new, i, a, k] <- get(algs[j])(dists, nk[k], hc.method)
          } else {
            arr[ind.new, i, a, k] <- get(algs[j])(dists, nk[k])
          }
        }
      }
    }
  }
  arr
}

#' Cluster other algorithms
#' @noRd
cc_other <- function(data, nk, p.item, reps, algs, xdim, ydim, rlen, alpha,
                     minPts, hc.method, seed.data, prep.data, scale, type,
                     min.var, pb, lalg, n) {
  alg <- toupper(algs)
  arr <- init_array(data, reps, alg, nk)

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
        arr[ind.new, i, j, k] <-
          switch(algs[j],
                 km = km(x, nk[k]),
                 ap = ap(x, nk[k]),
                 sc = sc(x, nk[k]),
                 gmm = gmm(x, nk[k]),
                 block = block(x, nk[k]),
                 som = som(x, nk[k], xdim, ydim, rlen, alpha, hc.method),
                 cmeans = cmeans(x, nk[k]),
                 hdbscan = hdbscan(x, minPts)
          )
      }
    }
  }
  arr
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
cdist <- function(x, dist, abs) {
  # Change partial matching distance methods from stats::dist to full names
  M <- c("euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")
  dist <- ifelse(dist %pin% M, M[pmatch(dist, M)], dist)

  # If spearman is specified, use spearman_dist()
  if (dist %pin% "spearman") return(spearman_dist(x = x, abs = abs))

  # If pearson is specified, use pearson_dist()
  if (dist %pin% "pearson") return(pearson_dist(x = x, abs = abs))

  # Identify custom distance functions from those found in stats::dist
  d <- try(stats::dist(x = x, method = dist), silent = TRUE)

  # Run custom function in parent environments if not found
  if (inherits(d, "try-error")) get(dist)(x) else d
}

#' Calculate pairwise Spearman correlational distances using
#' bioDist::spearman.dist defaults
#' @references
#' https://github.com/Bioconductor/bioDist/blob/devel/R/spearman.dist.R
#' @noRd
spearman_dist <- function(x, abs = TRUE) {
  if (abs) {
    scd <- 1 - abs(stats::cor(t(x), method = "spearman"))
  } else {
    scd <- 1 - stats::cor(t(x), method = "spearman")
  }
  scd %>%
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

#' Pearson distance
#' @noRd
pearson_dist <- function(x, abs = TRUE) {
  if (abs) {
    pd <- 1 - abs(stats::cor(t(x), method = "pearson"))
  } else {
    pd <- 1 - stats::cor(t(x), method = "pearson")
  }
  pd %>%
    magrittr::extract(lower.tri(.)) %>%
    `attributes<-`(
      list(
        Size = nrow(x),
        Labels = rownames(x),
        Diag = FALSE,
        Upper = FALSE,
        methods = "pearson",
        class = "dist"
      )
    )
}

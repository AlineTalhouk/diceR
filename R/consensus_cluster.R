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
  # Check for invalid distance inputs
  prep.data <- match.arg(prep.data)
  if (prep.data == "full")
    data <- prepare_data(data, scale = scale, type = type, min.var = min.var)

  # Use all algorithms if none are specified
  algorithms <- algorithms %||% ALG_NAMES

  # Store consensus dimensions for calculating progress bar increments/offsets
  lnk <- length(nk)
  lnmf <- sum(NALG %in% algorithms) * length(nmf.method)
  ldist <- sum(DALG %in% algorithms) * length(distance)
  lother <- sum(OALG %in% algorithms)
  if (progress) {
    pb <- utils::txtProgressBar(max = lnk * (lnmf + lother + ldist) * reps,
                                style = 3)
  } else {
    pb <- NULL
  }

  # Cluster NMF-based algorithms
  if (lnmf > 0) {
    nmf.arr <- cluster_nmf(data, nk, p.item, reps, nmf.method, seed.nmf,
                           seed.data, prep.data, scale, type, min.var,
                           progress, pb)
  } else {
    nmf.arr <- NULL
  }

  # Cluster distance-based algorithms
  dalgs <- algorithms[algorithms %in% DALG]
  if (ldist > 0) {
    dist.arr <- cluster_dist(data, nk, p.item, reps, dalgs, distance,
                             seed.data, prep.data, scale, type, min.var,
                             progress, pb, offset = lnk * lnmf * reps)
  } else {
    dist.arr <- NULL
  }

  # Cluster other algorithms
  oalgs <- algorithms[algorithms %in% OALG]
  if (lother > 0) {
    other.arr <- cluster_other(data, nk, p.item, reps, oalgs, xdim, ydim, rlen,
                               alpha, seed.data, prep.data, scale, type,
                               min.var, progress, pb, minPts,
                               offset = lnk * (lnmf + ldist) * reps)
    if ("hdbscan" %in% oalgs) {
      h.idx <- match("HDBSCAN", dimnames(other.arr)[[3]])
      h.obj <- other.arr[, , h.idx, ] %>%
        as.data.frame() %>%
        purrr::map(~ {
          c(prop_outlier = sum(.x == 0, na.rm = TRUE) / sum(!is.na(.x)),
            num_cluster = dplyr::n_distinct(!.x %in% c(NA, 0)))
        }) %>%
        purrr::transpose() %>%
        purrr::map(unlist)
      other.arr <- other.arr[, , -h.idx, , drop = FALSE]
    }
  } else {
    other.arr <- NULL
  }

  # Combine on third dimension (algorithm) and (optionally) save
  all.arr <- abind::abind(nmf.arr, dist.arr, other.arr, along = 3)
  if ("hdbscan" %in% algorithms) attr(all.arr, "hdbscan") <- h.obj
  if (!is.null(file.name)) {
    if (time.saved) {
      path <- paste0(file.name, "_",
                     format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), ".rds")
    } else {
      path <- paste0(file.name, ".rds")
    }
    saveRDS(all.arr, file = path)
  }
  all.arr
}

#' Cluster NMF-based algorithms
#' @noRd
cluster_nmf <- function(data, nk, p.item, reps, nmf.method, seed.nmf, seed.data,
                        prep.data, scale, type, min.var, progress, pb) {
  x_nmf <- nmf_transform(data)
  n <- nrow(data)
  lnmf <- length(nmf.method)
  lnk <- length(nk)
  nmf.arr <- array(NA, c(n, reps, lnmf, lnk),
                   dimnames = list(
                     rownames(data),
                     paste0("R", seq_len(reps)),
                     paste0("NMF_", Hmisc::capitalize(nmf.method)),
                     nk))
  for (k in seq_len(lnk)) {
    for (j in seq_len(lnmf)) {
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
        # Transpose since input for NMF::nmf uses rows as vars, cols as samples
        nmf.arr[ind.new, i, j, k] <- NMF::predict(NMF::nmf(
          t(x), rank = nk[k], method = nmf.method[j], seed = seed.nmf))
        if (progress) {
          utils::setTxtProgressBar(pb,
                                   (k - 1) * lnmf * reps + (j - 1) * reps + i)
        }
      }
    }
  }
  nmf.arr
}

#' Cluster algorithms with dissimilarity specification
#' @noRd
cluster_dist <- function(data, nk, p.item, reps, dalgs, distance, seed.data,
                         prep.data, scale, type, min.var, progress, pb,
                         offset) {
  n <- nrow(data)
  ld <- length(distance)
  lalg <- length(dalgs)
  ldist <- lalg * ld
  lnk <- length(nk)
  dist.arr <- array(NA, c(n, reps, ldist, lnk),
                    dimnames = list(
                      rownames(data),
                      paste0("R", seq_len(reps)),
                      apply(expand.grid(Hmisc::capitalize(distance),
                                        toupper(dalgs)),
                            1, function(x) paste0(x[2], "_", x[1])),
                      nk))
  for (k in seq_len(lnk)) {
    for (j in seq_len(lalg)) {
      for (d in seq_len(ld)) {
        set.seed(seed.data)
        for (i in seq_len(reps)) {
          # Find custom functions use get()
          ind.new <- sample(n, floor(n * p.item))
          if (prep.data == "sampled") {
            x <- prepare_data(data[ind.new, ], scale = scale, type = type,
                              min.var = min.var)
          } else if (prep.data %in% c("full", "none")) {
            x <- data[ind.new, ]
          }
          dists <- distances(x, distance[d])
          dist.arr[ind.new, i, (j - 1) * ld + d, k] <- get(dalgs[j])(dists[[1]],
                                                                     nk[k])
          if (progress)
            utils::setTxtProgressBar(pb, (k - 1) * lalg * ld * reps +
                                       (j - 1) * ld * reps +
                                       (d - 1) * reps + i + offset)
        }
      }
    }
  }
  dist.arr
}

#' Cluster other algorithms
#' @noRd
cluster_other <- function(data, nk, p.item, reps, oalgs, xdim, ydim, rlen,
                          alpha, seed.data, prep.data, scale, type, min.var,
                          progress, pb, minPts, offset) {
  n <- nrow(data)
  lalg <- length(oalgs)
  lnk <- length(nk)
  other.arr <- array(NA, c(n, reps, lalg, lnk),
                     dimnames = list(rownames(data),
                                     paste0("R", seq_len(reps)),
                                     toupper(oalgs),
                                     nk))
  for (k in seq_len(lnk)) {
    for (j in seq_len(lalg)) {
      set.seed(seed.data)
      for (i in seq_len(reps)) {
        ind.new <- sample(n, floor(n * p.item))
        if (prep.data == "sampled") {
          x <- prepare_data(data[ind.new, ], scale = scale, type = type,
                            min.var = min.var)
        } else if (prep.data %in% c("full", "none")) {
          x <- data[ind.new, ]
        }
        other.arr[ind.new, i, j, k] <-
          switch(oalgs[j],
                 ap = {
                   ap.cl <- stats::setNames(dplyr::dense_rank(suppressWarnings(
                     apcluster::apclusterK(apcluster::negDistMat, x,
                                           nk[k], verbose = FALSE)@idx)),
                     rownames(data[ind.new, ]))
                   if (length(ap.cl) == 0) NA else ap.cl
                 },
                 sc = stats::setNames(kernlab::specc(as.matrix(x), nk[k],
                                                     kernel = "rbfdot")@.Data,
                                      rownames(data[ind.new, ])),
                 gmm = mclust::Mclust(x, nk[k], verbose = FALSE)$classification,
                 block = {
                   blk.cl <- tryCatch(blockcluster::cocluster(
                     x, "continuous",
                     nbcocluster = c(nk[k], nk[k]))@rowclass + 1,
                     error = function(e) return(NA))
                   if (length(blk.cl) == 0) NA else blk.cl
                 },
                 som = som(x, nk[k], xdim = xdim, ydim = ydim, rlen = rlen,
                           alpha = alpha),
                 cmeans = cmeans(x, nk[k]),
                 hdbscan = dbscan::hdbscan(x = x, minPts = minPts)$cluster
          )
        if (progress)
          utils::setTxtProgressBar(pb, (k - 1) * lalg * reps +
                                     (j - 1) * reps + i + offset)
      }
    }
  }
  other.arr
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
    spear <- stats::setNames(list(spearman_dist(x)), "spearman")
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
  rvec <- stats::cor(t(x), method = "spearman") %>%
    abs() %>%
    magrittr::subtract(1, .) %>%
    magrittr::extract(lower.tri(.))
  attributes(rvec) <- list(Size = nrow(x), Labels = rownames(x), Diag = FALSE,
                           Upper = FALSE, methods = "spearman", class = "dist")
  rvec
}

#' Consensus clustering
#'
#' Runs consensus clustering across subsamples of the data, clustering
#' algorithms, and cluster sizes.
#'
#' See examples for how to use custom algorithms and distance functions. The
#' default clustering algorithms provided are:
#' \itemize{
#'   \item{"nmf": }{Nonnegative Matrix Factorization using Kullback-Leibler
#'   Divergence}
#'   \item{"nmfEucl": }{Nonnegative Matrix Factorization using Euclidean
#'   distance}
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
#' The \code{nmf.method} defaults are "brunet" (Kullback-Leibler divergence) and
#' "lee" (Euclidean distance).
#'
#' The progress bar increments for every unit of \code{reps}.
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
#' @param minPts number of minimum points in the eps region (for core points)
#'   for HDBSCAN. Default is 2 points.
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
#' @param save logical; if \code{TRUE}, the returned object will be saved at
#'   each iteration as well as at the end.
#' @param file.name file name of the written object
#' @param time.saved logical; if \code{TRUE}, the date saved is appended to the
#'   file name. Only applicable when \code{dir} is not \code{NULL}.
#' @return An array of dimension \code{nrow(x)} by \code{reps} by
#'   \code{length(algorithms)} by \code{length(nk)}. Each cube of the array
#'   represents a different k. Each slice of a cube is a matrix showing
#'   consensus clustering results for algorithms. The matrices have a row for
#'   each sample, and a column for each subsample. Each entry represents a class
#'   membership.
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
                              xdim = 10, ydim = 10, rlen = 200,
                              alpha = c(0.05, 0.01), minPts = 2,
                              distance = "euclidean",
                              prep.data = c("none", "full", "sampled"),
                              scale = TRUE, type = c("conventional", "robust"),
                              min.var = 1, progress = TRUE,
                              seed.nmf = 123456, seed.data = 1, save = FALSE,
                              file.name = "CCOutput", time.saved = FALSE) {
  # Check for invalid distance inputs
  prep.data <- match.arg(prep.data)
  if (prep.data == "full")
    data <- prepare_data(data, scale = scale, type = type, min.var = min.var)
  nmf.arr <- other.arr <- dist.arr <- NULL

  # Use all algorithms if none are specified
  algorithms <- algorithms %||% c("nmf", "hc", "diana", "km", "pam", "ap", "sc",
                                  "gmm", "block", "som", "cmeans", "hdbscan")

  # Store consensus dimensions for calculating progress bar increments/offsets
  lnk <- length(nk)
  lnmf <- ifelse("nmf" %in% algorithms, length(nmf.method), 0)
  ldist <- sum(!algorithms %in% c("nmf", "ap", "sc", "gmm", "block", "som",
                                  "cmeans", "hdbscan")) * length(distance)
  lother <- sum(c("ap", "sc", "gmm", "block", "som", "cmeans", "hdbscan") %in%
                  algorithms)
  if (progress) {
    pb <- utils::txtProgressBar(max = lnk * (lnmf + lother + ldist) * reps,
                                style = 3)
  } else {
    pb <- NULL
  }

  # Cluster NMF-based algorithms
  if ("nmf" %in% algorithms) {
    nmf.arr <- cluster_nmf(data, nk, p.item, reps, nmf.method,
                           seed.nmf, seed.data, prep.data, scale, type, min.var,
                           progress, pb)
  }

  # Cluster distance-based algorithms
  dalgs <- algorithms[!algorithms %in% c("nmf", "ap", "sc", "gmm", "block",
                                         "som", "cmeans", "hdbscan")]
  if (length(dalgs) > 0) {
    dist.arr <- cluster_dist(data, nk, p.item, reps, dalgs, distance,
                             seed.data, prep.data, scale, type, min.var,
                             progress, pb, offset = lnk * lnmf * reps)
  }

  # Cluster other algorithms
  oalgs <- algorithms[algorithms %in% c("ap", "sc", "gmm", "block",
                                        "som", "cmeans", "hdbscan")]
  if (length(oalgs) > 0) {
    other.arr <- cluster_other(data, nk, p.item, reps, oalgs, xdim, ydim, rlen,
                               alpha, seed.data, prep.data, scale, type,
                               min.var, progress, pb, eps, minPts,
                               offset = lnk * (lnmf + ldist) * reps)
  }

  # Combine on third dimension (algorithm) and (optionally) save
  all.arr <- abind::abind(nmf.arr, dist.arr, other.arr, along = 3)
  if (save) {
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
  # Transform to non-negative matrix by column-binding a negative replicate and
  # then coercing all negatives to 0
  x.nmf <- data %>%
    cbind(-.) %>%
    apply(2, function(d) ifelse(d < 0, 0, d))
  n <- nrow(data)
  n.new <- floor(n * p.item)
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
        ind.new <- sample(n, n.new, replace = FALSE)
        # Transpose since input for NMF::nmf uses rows as vars, cols as samples
        # In case the subsample has all-zero vars, remove them to speed up comp
        x.nmf_samp <- x.nmf[ind.new, !(apply(x.nmf[ind.new, ], 2,
                                             function(x) all(x == 0)))]
        if (prep.data == "sampled") {
          x <- prepare_data(x.nmf_samp, scale = scale, type = type,
                            min.var = min.var) %>%
            cbind(-.) %>%
            apply(2, function(d) ifelse(d < 0, 0, d))
        } else if (prep.data %in% c("full", "none")) {
          x <- x.nmf_samp
        }
        nmf.arr[ind.new, i, j, k] <- NMF::predict(NMF::nmf(
          t(x), rank = nk[k], method = nmf.method[j], seed = seed.nmf))
        if (progress)
          utils::setTxtProgressBar(pb,
                                   (k - 1) * lnmf * reps + (j - 1) * reps + i)
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
  n.new <- floor(n * p.item)
  ld <- length(distance)
  lalg <- length(dalgs)
  ldist <- prod(lalg, ld)
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
          ind.new <- sample(n, n.new, replace = FALSE)
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
                          progress, pb, eps, minPts, offset) {
  n <- nrow(data)
  n.new <- floor(n * p.item)
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
        ind.new <- sample(n, n.new, replace = FALSE)
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

#' @noRd
hc <- function(d, k, method = "average") {
  as.integer(stats::cutree(stats::hclust(d, method = method), k))
}

#' @noRd
diana <- function(d, k) {
  as.integer(stats::cutree(cluster::diana(d, diss = TRUE), k))
}

#' @noRd
km <- function(d, k) {
  as.integer(stats::kmeans(d, k)$cluster)
}

#' @noRd
pam <- function(d, k) {
  as.integer(cluster::pam(d, k, cluster.only = TRUE))
}

#' Main function to train SOM and then do hierarchical clustering
#' @noRd
som <- function(x, k, xdim, ydim, rlen, alpha, method = "average") {
  if (!is.matrix(x)) x <- as.matrix(x)
  model <- som_train(x = x, xdim = xdim, ydim = ydim, rlen = rlen,
                     alpha = alpha)
  som_cluster(model = model, k = k, method = method)
}

#' Train the SOM, specifiy grid size and other optional parameters.
#' @noRd
som_train <- function(x, xdim, ydim, rlen, alpha, topo = "hexagonal") {
  # Create SOM grid and map data into the grid
  grid <- kohonen::somgrid(xdim = xdim, ydim = ydim, topo = topo)
  kohonen::som(x, grid = grid, rlen = rlen, alpha = alpha)
}

#' Cut tree into k groups and output cluster labels for original data
#' @noRd
som_cluster <- function(model, k, method) {
  # Get distance matrix, use hc to cluster the codebook vectors
  cl <- hc(stats::dist(kohonen::getCodes(model, 1)),
           k = k,
           method = method)
  pred <- stats::predict(model)$unit.classif
  cl[pred]
}

#' Fuzzy c-means using best m via validity/performance measures
#' @noRd
cmeans <- function(x, k) {
  fuzzy <- seq(1.1, 3, by = 0.1) %>%
    purrr::map(~ e1071::cmeans(x = x, centers = k, m = .x))
  mbest <- c("xie.beni", "fukuyama.sugeno", "partition.entropy") %>%
    purrr::map(
      ~ purrr::map(fuzzy,
                   function(f) e1071::fclustIndex(y = f, x = x, index = .x)
      )
    ) %>%
    purrr::map_int(which.min) %>%
    table() %>%
    which.max() %>%
    names() %>%
    as.integer()
  fuzzy[[mbest]]$cluster
}

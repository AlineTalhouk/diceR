#' Consensus clustering over samples and algorithms
#' 
#' Generates multiple runs for consensus clustering among replicated subsamples 
#' of a dataset as well as across different clustering algorithms.
#' 
#' The clustering algorithms provided are:
#' \itemize{
#'   \item{"nmf": }{Nonnegative Matrix Factorization using Kullback-Leibler Divergence} 
#'   \item{"nmfEucl": }{Nonnegative Matrix Factorization using Euclidean distance}
#'   \item{"hc": }{Hierarchical Clustering}
#'   \item{"diana": }{DIvisive ANAlysis Clustering}
#'   \item{"km": }{K-Means Clustering}
#'   \item{"pam": }{Partition Around Mediods}
#'   \item{"ap": }{Affinity Propagation}
#'   \item{"sc": }{Spectral Clustering using Radial-Basis kernel function}
#'   \item{"gmm": }{Gaussian Mixture Model using Bayesian Information Criterion on EM algorithm} 
#'   \item{"block": }{Biclustering using a latent block model}
#' }
#' 
#' The \code{nmf.method} defaults are "brunet" (Kullback-Leibler divergence) and
#' "lee" (Euclidean distance).
#' 
#' The \code{min.sd} argument is used to filter the feature space for only 
#' highly variable features. Only features with a standard deviation across all 
#' samples greater than \code{min.sd} will be used.
#' 
#' The progress bar increments for every unit of \code{reps}.
#' 
#' @param data data matrix with rows as samples and columns as variables
#' @param nk number of clusters (k) requested; can specify a single integer or a
#'   range of integers to compute multiple k
#' @param pItem proportion of items to be used in subsampling within an 
#'   algorithm
#' @param reps number of subsamples
#' @param algorithms vector of clustering algorithms for performing consensus 
#'   clustering. Must be any number of the following: "nmf", 
#'   "hc", "diana", "km", "pam", "ap", "sc", "gmm", "block". See details. Can
#'   use a custom clustering algorithm. See example.
#' @param nmf.method specify NMF-based algorithms to run. By default the 
#'   "brunet" and "lee" algorithms are called. See \code{\link[NMF]{nmf}} for
#'   details.
#' @param distance a vector of distance functions. Defaults to "euclidean". Can 
#'   use a custom distance function. See example.
#' @param progress logical; should a progress bar be displayed?
#' @param seed random seed to use for NMF-based algorithms
#' @param seed.alg seed to use to ensure each algorithm operates on the same set
#'   of subsamples
#' @param min.sd minimum standard deviation threshold. See details.
#' @param save logical; if \code{TRUE}, the returned object will be saved at 
#'   each iteration as well as at the end.
#' @param file.name file name of the written object
#' @param time.saved logical; if \code{TRUE}, the date saved is appended to the 
#'   file name. Only applicable when \code{dir} is not \code{NULL}.
#' @return An array of dimension \code{nrow(x)} by \code{reps} by 
#'   \code{length(algorithms)} Each slice of the array is a matrix showing 
#'   consensus clustering results for algorithms. The matrices have a row for
#'   each sample, and a column for each subsample. Each entry represents a class
#'   membership.
#' @author Derek Chiu, Aline Talhouk
#' @import mclust
#' @export
#' @examples 
#' data(hgsc)
#' dat <- t(hgsc[, -1])
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
#' cc <- consensus_cluster(dat, reps = 5, algorithms = c("pam", "agnes"),
#' distance = c("euclidean", "manh"))
#' str(cc)
consensus_cluster <- function(data, nk = 2:4, pItem = 0.8, reps = 1000,
                              algorithms = NULL, nmf.method = c("brunet", "lee"),
                              distance = "euclidean",
                              progress = TRUE, seed = 123456, seed.alg = 1,
                              min.sd = 1, save = FALSE, file.name = "CCOutput",
                              time.saved = FALSE) {
  data.prep <- prepare_data(data, min.sd = min.sd)
  nmf.arr <- other.arr <- dist.arr <- NULL
  check.dists <- distances(data.prep, distance)
  lnk <- length(nk)
  lnmf <- ifelse("nmf" %in% algorithms, length(nmf.method), 0)
  ldist <- sum(!algorithms %in% c("nmf", "ap", "sc", "gmm", "block")) * length(distance)
  lother <- sum(c("ap", "sc", "gmm", "block") %in% algorithms)
  if (progress) {
    pb <- utils::txtProgressBar(max = lnk * (lnmf + lother + ldist) * reps, style = 3)
  } else {
    pb <- NULL
  }
  
  if (is.null(algorithms))
    algorithms <- c("nmf", "hc", "diana", "km", "pam",
                    "ap", "sc", "gmm", "block")
  
  if ("nmf" %in% algorithms) {
    nmf.arr <- cluster_nmf(data.prep, nk, pItem, reps, nmf.method,
                           seed, seed.alg, progress, pb)
  }
  
  dalgs <- algorithms[!algorithms %in% c("nmf", "ap", "sc", "gmm", "block")]
  if (length(dalgs) > 0) {
    dist.arr <- cluster_dist(data.prep, nk, pItem, reps, dalgs, distance,
                             seed, seed.alg, progress, pb,
                             offset = lnk * lnmf * reps)
  }
  
  oalgs <- algorithms[algorithms %in% c("ap", "sc", "gmm", "block")]
  if (length(oalgs) > 0) {
    other.arr <- cluster_other(data.prep, nk, pItem, reps, oalgs,
                               seed, seed.alg, progress, pb,
                               offset = lnk * (lnmf + ldist) * reps)
  }
  
  all.arr <- abind::abind(nmf.arr, dist.arr, other.arr, along = 3)
  if (save) {
    if (time.saved) {
      path <- paste0(file.name, "_",
                     format(Sys.time(), "%Y-%m-%d_%H-%M-%S") , ".rds")
    } else {
      path <- paste0(file.name, ".rds")
    }
    readr::write_rds(all.arr, path = path)
  }
  return(all.arr)
}

#' Prepare data for consensus clustering
#'
#' Remove variables with low signal and scale before consensus clustering
#'
#' The \code{min.sd} argument is used to filter the feature space for only
#' highly variable features. Only features with a standard deviation across all
#' samples greater than \code{min.sd} will be used.
#'
#' @param data data matrix with rows as samples and columns as variables
#' @param min.sd minimum standard deviation threshold. See details.
#' @return dataset prepared for usage in \code{consensus_cluster}
#' @author Derek Chiu
#' @export
#' @examples
#' set.seed(2)
#' x <- replicate(10, rnorm(100))
#' x.prep <- prepare_data(x)
#' dim(x)
#' dim(x.prep)
prepare_data <- function(data, min.sd = 1) {
  dat.out <- data %>%
    magrittr::extract(apply(., 1, function(x) !any(is.na(x))),
                      apply(., 2, function(x) stats::sd(x, na.rm = TRUE)) >
                        min.sd) %>%
    scale()
  return(dat.out)
}

#' Cluster NMF-based algorithms
#' @noRd
cluster_nmf <- function(data, nk, pItem, reps, nmf.method, seed, seed.alg,
                        progress, pb) {
  x.nmf <- data %>%
    cbind(-.) %>%
    apply(2, function(x) ifelse(x < 0, 0, x))
  n <- nrow(data)
  n.new <- floor(n * pItem)
  lnmf <- length(nmf.method)
  lnk <- length(nk)
  nmf.arr <- array(NA, c(n, reps, lnmf, lnk),
                   dimnames = list(rownames(data),
                                   paste0("R", 1:reps),
                                   paste0("NMF_", Hmisc::capitalize(nmf.method)),
                                   nk))
  for (k in 1:lnk) {
    for (j in 1:lnmf) {
      set.seed(seed.alg)
      for (i in 1:reps) {
        ind.new <- sample(n, n.new, replace = FALSE)
        x.nmf.samp <- t(x.nmf[ind.new, !(apply(x.nmf[ind.new, ], 2,
                                               function(x) all(x == 0)))])
        nmf.arr[ind.new, i, j, k] <- NMF::predict(NMF::nmf(
          x.nmf.samp, rank = nk[k], method = nmf.method[j], seed = seed))
        if (progress)
          utils::setTxtProgressBar(pb, (k - 1) * lnmf * reps + (j - 1) * reps + i)
      }
    }
  }
  return(nmf.arr)
}

#' Cluster algorithms with dissimilarity specification
#' @noRd
cluster_dist <- function(data, nk, pItem, reps, dalgs, distance, seed, seed.alg,
                         progress, pb, offset) {
  n <- nrow(data)
  n.new <- floor(n * pItem)
  ld <- length(distance)
  lalg <- length(dalgs)
  ldist <- prod(lalg, ld)
  lnk <- length(nk)
  dist.arr <- array(NA, c(n, reps, ldist, lnk),
                    dimnames = list(rownames(data),
                                    paste0("R", 1:reps),
                                    apply(expand.grid(Hmisc::capitalize(distance),
                                                      toupper(dalgs)),
                                          1, function(x) paste0(x[2], "_", x[1])),
                                    nk))
  for (k in 1:lnk) {
    for (j in 1:lalg) {
      for (d in 1:ld) {
        set.seed(seed.alg)
        for (i in 1:reps) {
          ind.new <- sample(n, n.new, replace = FALSE)
          dists <- distances(data[ind.new, ], distance[d])
          dist.arr[ind.new, i, (j - 1) * ld + d, k] <- get(dalgs[j])(dists[[1]], nk[k])
          if (progress)
            utils::setTxtProgressBar(pb, (k - 1) * lalg * ld * reps +
                                       (j - 1) * ld * reps +
                                       (d - 1) * reps + i + offset)
        }
      }
    }
  }
  return(dist.arr)
}

#' Cluster other algorithms
#' @noRd
cluster_other <- function(data, nk, pItem, reps, oalgs, seed, seed.alg,
                          progress, pb, offset) {
  n <- nrow(data)
  n.new <- floor(n * pItem)
  lalg <- length(oalgs)
  lnk <- length(nk)
  other.arr <- array(NA, c(n, reps, lalg, lnk),
                     dimnames = list(rownames(data),
                                     paste0("R", 1:reps),
                                     toupper(oalgs),
                                     nk))
  for (k in 1:lnk) {
    for (j in 1:lalg) {
      set.seed(seed.alg)
      for (i in 1:reps) {
        ind.new <- sample(n, n.new, replace = FALSE)
        other.arr[ind.new, i, j, k] <- 
          switch(oalgs[j],
                 ap = stats::setNames(dplyr::dense_rank(suppressWarnings(
                   apcluster::apclusterK(apcluster::negDistMat, data[ind.new, ],
                                         nk[k], verbose = FALSE)@idx)),
                   rownames(data[ind.new, ])),
                 sc = stats::setNames(kernlab::specc(data[ind.new, ], nk[k],
                                                     kernel = "rbfdot")@.Data,
                                      rownames(data[ind.new, ])),
                 gmm = mclust::Mclust(data[ind.new, ], nk[k])$classification,
                 block = blockcluster::cocluster(
                   data[ind.new, ], "continuous",
                   nbcocluster = c(nk[k], nk[k]))@rowclass + 1)
        if (progress)
          utils::setTxtProgressBar(pb, (k - 1) * lalg * reps +
                                     (j - 1) * reps + i + offset)
      }
    }
  }
  return(other.arr)
}

#' Return a list of distance matrices
#' @param x data matrix
#' @param dist a character vector of distance methods taken from stats::dist, or
#'   "spearman", or a custom distance function in the current environment
#' @noRd
distances <- function(x, dist) {
  # Change partial matching distance methods from stats::dist to full names
  METHODS <- c("euclidean", "maximum", "manhattan", "canberra",
               "binary", "minkowski")
  dist <- ifelse(dist %pin% METHODS, METHODS[pmatch(dist, METHODS)], dist)
  
  # Check if spearman distance is used (starts with string)
  sp.idx <- sapply(paste0("^", dist), grepl, "spearman")
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
    check <- stats::setNames(lapply(d, function(d)
      try(stats::dist(x = x, method = d), silent = TRUE)), d)
    is.error <- sapply(check, inherits, "try-error")
    succeeded <- which(!is.error)
    failed <- which(is.error)
    
    # Search for custom function in parent environments
    if (length(failed) > 0) {
      custom <- stats::setNames(lapply(names(check[failed]),
                                       function(d) get(d)(x)),
                                names(failed))
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
#' @references https://github.com/Bioconductor-mirror/bioDist/blob/master/R/spearman.dist.R
#' @noRd
spearman_dist <- function(x) {
  rvec <- stats::cor(t(x), method = "spearman") %>% 
    abs() %>% 
    magrittr::subtract(1, .) %>% 
    magrittr::extract(lower.tri(.))
  attributes(rvec) <- list(Size = nrow(x), Labels = rownames(x), Diag = FALSE,
                           Upper = FALSE, methods = "spearman", class = "dist")
  return(rvec)
}

#' @noRd
hc <- function(d, k, method = "average") {
  return(as.integer(stats::cutree(stats::hclust(d, method = method), k)))
}

#' @noRd
diana <- function(d, k) {
  return(as.integer(stats::cutree(cluster::diana(d, diss = TRUE), k)))
}

#' @noRd
km <- function(d, k) {
  return(as.integer(stats::kmeans(d, k)$cluster))
}

#' @noRd
pam <- function(d, k) {
  return(as.integer(cluster::pam(d, k, cluster.only = TRUE)))
}
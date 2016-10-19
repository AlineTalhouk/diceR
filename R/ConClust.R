#' Consensus clustering over samples and algorithms
#' 
#' Generates multiple runs for consensus clustering among replicated subsamples 
#' of a dataset as well as across different clustering algorithms.
#' 
#' The clustering algorithms provided in \code{ConClust} are:
#' \itemize{ 
#' \item{"nmfDiv": }{Nonnegative Matrix Factorization using Kullback-Leibler Divergence}
#' \item{"nmfEucl": }{Nonnegative Matrix Factorization using Euclidean distance}
#' \item{"hcAEucl": }{Hierarchical Clustering using Average linkage and Euclidean distance}
#' \item{"hcDianaEucl": }{Divisive Hierarchical Clustering using Euclidean distance}
#' \item{"kmEucl": }{K-Means using Euclidean distance}
#' \item{"kmSpear": }{K-Means using Spearman distance} 
#' \item{"pamEucl": }{Partition Around Mediods using Euclidean distance} 
#' \item{"pamSpear": }{Partition Around Mediods using Spearman distance} 
#' \item{"apEucl": }{Affinity Propagation using Euclidean distance} 
#' \item{"scRbf": }{Spectral Clustering using Radial-Basis kernel function} 
#' \item{"gmmBIC": }{Gaussian Mixture Model using Bayesian Information Criterion on EM algorithm}
#' \item{"biclust": }{Biclustering using a latent block model} 
#' }
#' 
#' The \code{min.sd} argument is used to filter the feature space for only 
#' highly variable features. Only features with a standard deviation across all 
#' samples greater than \code{min.sd} will be used.
#' 
#' The progress bar for \code{parallel = FALSE} increments for every unit of 
#' \code{reps}, where as for \code{parallel = TRUE}, the increments are for 
#' every \code{method}. Additionally, when the progress bar is used in 
#' conjunction with parallel computation, each worker outputs the attached 
#' packages (and potential masking of functions).
#' 
#' @param x data matrix with rows as samples and columns as variables
#' @param k number of clusters requested
#' @param pItem proportion of items to be used in subsampling within an 
#'   algorithm
#' @param reps number of subsamples
#' @param method vector of clustering algorithms for performing consensus 
#'   clustering. Must be any number of the following: "nmfDiv", "nmfEucl", 
#'   "hcAEucl", "hcDianaEucl", "kmEucl", "kmSpear", "pamEucl", "pamSpear", 
#'   "apEucl", "scRbf", "gmmBIC", "biclust". See details.
#' @param parallel logical; if \code{TRUE}, the function registers a parallel 
#'   backend from the \code{doParallel} package and runs a foreach loop. By 
#'   default, \code{parallel = NULL} means the function determines, based on 
#'   other parameters, whether to use parallel or not (e.g. number of
#'   \code{reps}, availability of multiple cores).
#' @param ncores number of CPU cores to use for parallel computation. Default 
#'   uses all available cores.
#' @param progress logical; should a progress bar be displayed?
#' @param seed random seed to use for NMF-based algorithms
#' @param seed.method seed to use to ensure each method operates on the same set
#'   of subsamples
#' @param min.sd minimum standard deviation threshold. See details.
#' @param save logical; if \code{TRUE}, the returned object will be saved at 
#'   each iteration as well as at the end.
#' @param file.name file name of the written object
#' @param time.saved logical; if \code{TRUE}, the date saved is appended to the 
#'   file name. Only applicable when \code{dir} is not \code{NULL}.
#' @return An array of dimension \code{nrow(x)} by \code{reps} by 
#'   \code{length(methods)} Each slice of the array is a matrix showing 
#'   consensus clustering results for algorithms (method). The matrices have a 
#'   row for each sample, and a column for each subsample. Each entry represents
#'   a class membership.
#' @author Derek Chiu, Aline Talhouk
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats dist hclust cutree kmeans setNames
#' @import foreach mclust
#' @export
#' @examples
#' data(hgsc)
#' dat <- t(hgsc[, -1])
#' x <- ConClust(dat, k = 4, reps = 10, method = "hcAEucl", save = FALSE)
ConClust <- function(x, k, pItem = 0.8, reps = 1000, method = NULL,
                     parallel = NULL, ncores = NULL, progress = TRUE,
                     seed = 123456, seed.method = 1, min.sd = 1, save = FALSE,
                     file.name = "ConClustOutput", time.saved = FALSE) {
  if (is.null(method))
    method <- c("nmfDiv", "nmfEucl", "hcAEucl", "hcDianaEucl", "kmEucl",
                "kmSpear", "pamEucl", "pamSpear", "apEucl", "scRbf", "gmmBIC",
                "biclust")
  if (is.null(ncores)) ncores <- parallel::detectCores()
  if (is.null(parallel)) {
    if (reps >= 100 & ncores >= 2) {
      dopar <- TRUE
    } else {
      dopar <- FALSE
    }
  } else if (isTRUE(parallel)) {
    dopar <- TRUE
  } else if (!isTRUE(parallel)) {
    dopar <- FALSE
  }
  x.rest <- prepare_data(x, min.sd = min.sd)
  x.nmf <- x.rest %>%
    cbind(-.) %>%
    apply(2, function(x) ifelse(x < 0, 0, x))
  samples <- rownames(x.rest)
  n <- nrow(x.rest)
  n.new <- floor(n * pItem)
  nm <- length(method)
  coclus <- array(NA, c(n, reps, nm),
                  dimnames = list(samples, paste0("R", 1:reps), method))
  if (!dopar) {
    if (progress) pb <- txtProgressBar(max = nm * reps, style = 3)
    for (j in 1:nm) {
      set.seed(seed.method)
      for (i in 1:reps) {
        if (progress) setTxtProgressBar(pb, (j - 1) * reps + i)
        ind.new <- sample(n, n.new, replace = FALSE)
        if (any(c("nmfDiv", "nmfEucl") %in% method))
          x.nmf.samp <- x.nmf[ind.new, !(apply(x.nmf[ind.new, ], 2,
                                               function(x) all(x == 0)))]
        coclus[ind.new, i, j] <- switch(
          method[j],
          nmfDiv = NMF::predict(NMF::nmf(
            t(x.nmf.samp), rank = k, method = "brunet", seed = seed)),
          nmfEucl = NMF::predict(NMF::nmf(
            t(x.nmf.samp), rank = k, method = "lee", seed = seed)),
          hcAEucl = cutree(hclust(dist(x.rest[ind.new, ]),
                                  method = "average"), k),
          hcDianaEucl = cutree(cluster::diana(euc(x.rest[ind.new, ]),
                                              diss = TRUE), k),
          kmEucl = kmeans(euc(x.rest[ind.new, ]),
                          k)$cluster,
          kmSpear = kmeans(spearman.dist(x.rest[ind.new, ]),
                           k)$cluster,
          pamEucl = cluster::pam(euc(x.rest[ind.new, ]), k,
                                 cluster.only = TRUE),
          pamSpear = cluster::pam(spearman.dist(x.rest[ind.new, ]), k,
                                  cluster.only = TRUE),
          apEucl = setNames(dense_rank(suppressWarnings(
            apcluster::apclusterK(apcluster::negDistMat, x.rest[ind.new, ], k,
                                  verbose = FALSE)@idx)),
            rownames(x.rest[ind.new, ])),
          scRbf = setNames(kernlab::specc(x.rest[ind.new, ], k,
                                          kernel = "rbfdot")@.Data,
                           rownames(x.rest[ind.new, ])),
          gmmBIC = mclust::Mclust(x.rest[ind.new, ], k)$classification,
          biclust = blockcluster::cocluster(x.rest[ind.new, ],
                                            "continuous",
                                            nbcocluster = c(k, k))@rowclass + 1
        )
        if (i %% 10 == 0 & save) {
          if (time.saved) {
            path <- paste0(file.name, "_",
                           format(Sys.time(), "%Y-%m-%d_%H-%M-%S") , ".rds")
          } else {
            path <- paste0(file.name, ".rds")
          }
          readr::write_rds(coclus, path = path)
        }
      }
    }
  } else {
    if (progress) {
      cl <- parallel::makeCluster(ncores, outfile = "")
    } else {
      cl <- parallel::makeCluster(ncores)
    }
    doParallel::registerDoParallel(cl)
    coclus <- foreach(j = 1:nm, .packages = c("foreach", "bioDist", "dplyr", "apcluster", "mclust")) %dopar% {
      set.seed(seed.method)
      if (progress) pb <- txtProgressBar(max = reps, style = 3)
      foreach(i = 1:reps) %dopar% {
        if (progress) setTxtProgressBar(pb, (j - 1) * reps + i)
        ind.new <- sample(n, n.new, replace = FALSE)
        if (any(c("nmfDiv", "nmfEucl") %in% method))
          x.nmf.samp <- x.nmf[ind.new, !(apply(x.nmf[ind.new, ], 2,
                                               function(x) all(x == 0)))]
        coclus[ind.new, i, j] <- switch(
          method[j],
          nmfDiv = NMF::predict(NMF::nmf(
            t(x.nmf.samp), rank = k, method = "brunet", seed = seed)),
          nmfEucl = NMF::predict(NMF::nmf(
            t(x.nmf.samp), rank = k, method = "lee", seed = seed)),
          hcAEucl = cutree(hclust(dist(x.rest[ind.new, ]),
                                  method = "average"), k),
          hcDianaEucl = cutree(cluster::diana(euc(x.rest[ind.new, ]),
                                              diss = TRUE), k),
          kmEucl = kmeans(euc(x.rest[ind.new, ]),
                          k)$cluster,
          kmSpear = kmeans(spearman.dist(x.rest[ind.new, ]),
                           k)$cluster,
          pamEucl = cluster::pam(euc(x.rest[ind.new, ]), k,
                                 cluster.only = TRUE),
          pamSpear = cluster::pam(spearman.dist(x.rest[ind.new, ]), k,
                                  cluster.only = TRUE),
          apEucl = setNames(dense_rank(suppressWarnings(
            apcluster::apclusterK(apcluster::negDistMat, x.rest[ind.new, ], k,
                                  verbose = FALSE)@idx)),
            rownames(x.rest[ind.new, ])),
          scRbf = setNames(kernlab::specc(x.rest[ind.new, ], k,
                                          kernel = "rbfdot")@.Data,
                           rownames(x.rest[ind.new, ])),
          gmmBIC = mclust::Mclust(x.rest[ind.new, ], k)$classification,
          biclust = blockcluster::cocluster(x.rest[ind.new, ],
                                            "continuous",
                                            nbcocluster = c(k, k))@rowclass + 1
        )
      }
      coclus[, , j]
    } %>%
      unlist() %>% 
      array(dim = c(n, reps, nm),
            dimnames = list(samples, paste0("R", 1:reps), method))
    doParallel::stopImplicitCluster()
  }
  if (save) {
    if (time.saved) {
      path <- paste0(file.name, "_",
                     format(Sys.time(), "%Y-%m-%d_%H-%M-%S") , ".rds")
    } else {
      path <- paste0(file.name, ".rds")
    }
    readr::write_rds(coclus, path = path)
  }
  return(coclus)
}

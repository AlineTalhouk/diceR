#' Consensus clustering over samples and algorithms
#' 
#' Generates multiple runs for consensus clustering among replicated subsamples 
#' of a dataset as well as across different clustering algorithms.
#' 
#' The clustering algorithms provided are:
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
#' \code{algorithms}. Additionally, when the progress bar is used in 
#' conjunction with parallel computation, each worker outputs the attached 
#' packages (and potential masking of functions).
#' 
#' @param data data matrix with rows as samples and columns as variables
#' @param nk number of clusters (k) requested; can specify a single integer or a
#'   range of integers to compute multiple k
#' @param pItem proportion of items to be used in subsampling within an 
#'   algorithm
#' @param reps number of subsamples
#' @param algorithms vector of clustering algorithms for performing consensus 
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
#'   consensus clustering results for algorithms. The matrices have a 
#'   row for each sample, and a column for each subsample. Each entry represents
#'   a class membership.
#' @author Derek Chiu, Aline Talhouk
#' @import foreach mclust
#' @export
#' @examples
#' data(hgsc)
#' dat <- t(hgsc[, -1])
#' x1 <- consensus_cluster(dat, nk = 2:4, reps = 10, algorithms = c("hcAEucl",
#' "kmEucl"))
consensus_cluster <- function(data, nk = 2:4, pItem = 0.8, reps = 1000,
                              algorithms = NULL, parallel = NULL, ncores = NULL,
                              progress = TRUE, seed = 123456, seed.alg = 1,
                              min.sd = 1, save = FALSE, file.name = "CCOutput",
                              time.saved = FALSE) {
  if (is.null(algorithms))
    algorithms <- c("nmfDiv", "nmfEucl", "hcAEucl", "hcDianaEucl", "kmEucl",
                    "kmSpear", "pamEucl", "pamSpear", "apEucl", "scRbf",
                    "gmmBIC", "biclust")
  if (is.null(ncores)) ncores <- parallel::detectCores() - 1
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
  x.rest <- prepare_data(data, min.sd = min.sd)
  x.nmf <- x.rest %>%
    cbind(-.) %>%
    apply(2, function(x) ifelse(x < 0, 0, x))
  samples <- rownames(x.rest)
  n <- nrow(x.rest)
  n.new <- floor(n * pItem)
  lnk <- length(nk)
  nm <- length(algorithms)
  coclus <- array(NA, c(n, reps, nm, lnk),
                  dimnames = list(samples, paste0("R", 1:reps), algorithms, nk))
  if (!dopar) {
    if (progress)
      pb <- utils::txtProgressBar(max = lnk * nm * reps, style = 3)
    for (k in 1:lnk) {
      for (j in 1:nm) {
        set.seed(seed.alg)
        for (i in 1:reps) {
          if (progress) 
            utils::setTxtProgressBar(pb,
                                     (k - 1) * nm * reps + (j - 1) * reps + i)
          ind.new <- sample(n, n.new, replace = FALSE)
          if (any(c("nmfDiv", "nmfEucl") %in% algorithms))
            x.nmf.samp <- t(x.nmf[ind.new, !(apply(x.nmf[ind.new, ], 2,
                                                 function(x) all(x == 0)))])
          coclus[ind.new, i, j, k] <- switch(
            algorithms[j],
            nmfDiv = NMF::predict(NMF::nmf(
              x.nmf.samp, rank = nk[k], method = "brunet", seed = seed)),
            nmfEucl = NMF::predict(NMF::nmf(
              x.nmf.samp, rank = nk[k], method = "lee", seed = seed)),
            hcAEucl = hcAEucl(x.rest[ind.new, ], nk[k]),
            hcDianaEucl = hcDianaEucl(x.rest[ind.new, ], nk[k]),
            kmEucl = stats::kmeans(
              stats::dist(x.rest[ind.new, ]), nk[k])$cluster,
            kmSpear = stats::kmeans(
              spearman_dist(x.rest[ind.new, ]), nk[k])$cluster,
            pamEucl = cluster::pam(
              stats::dist(x.rest[ind.new, ]), nk[k], cluster.only = TRUE),
            pamSpear = cluster::pam(
              spearman_dist(x.rest[ind.new, ]), nk[k], cluster.only = TRUE),
            apEucl = stats::setNames(dplyr::dense_rank(suppressWarnings(
              apcluster::apclusterK(apcluster::negDistMat, x.rest[ind.new, ],
                                    nk[k], verbose = FALSE)@idx)),
              rownames(x.rest[ind.new, ])),
            scRbf = stats::setNames(kernlab::specc(x.rest[ind.new, ], nk[k],
                                                   kernel = "rbfdot")@.Data,
                                    rownames(x.rest[ind.new, ])),
            gmmBIC = mclust::Mclust(x.rest[ind.new, ], nk[k])$classification,
            biclust = blockcluster::cocluster(
              x.rest[ind.new, ], "continuous",
              nbcocluster = c(nk[k], nk[k]))@rowclass + 1
          )
          if (i %% 10 == 0 & save) {
            if (time.saved) {
              path <- paste0(file.name, "_",
                             format(Sys.time(), "%Y-%m-%d_%H-%M-%S") , ".rds")
            } else {
              path <- paste0(file.name, ".rds")
            }
            readr::write_rds(coclus, path = path)
            message(paste("Second dimension of coclus is:",
                          dim(readRDS(path))[2]))
          }
        }
      }
    }
  } else {
    if (progress) {
      cl <- parallel::makeCluster(ncores, useXDR = FALSE, outfile = "")
    } else {
      cl <- parallel::makeCluster(ncores, useXDR = FALSE)
    }
    doParallel::registerDoParallel(cl)
    coclus <- foreach(k = 1:lnk, .packages = c("foreach", "bioDist", "dplyr", "apcluster", "mclust")) %dopar% {
      foreach(j = 1:nm) %dopar% {
        set.seed(seed.alg)
        if (progress)
          pb <- utils::txtProgressBar(max = nm * reps, style = 3)
        foreach(i = 1:reps) %dopar% {
          if (progress)
            utils::setTxtProgressBar(pb,
                                     (k - 1) * nm * reps + (j - 1) * reps + i)
          ind.new <- sample(n, n.new, replace = FALSE)
          if (any(c("nmfDiv", "nmfEucl") %in% algorithms))
            x.nmf.samp <- t(x.nmf[ind.new, !(apply(x.nmf[ind.new, ], 2,
                                                 function(x) all(x == 0)))])
          coclus[ind.new, i, j, k] <- switch(
            algorithms[j],
            nmfDiv = NMF::predict(NMF::nmf(
              x.nmf.samp, rank = nk[k], method = "brunet", seed = seed)),
            nmfEucl = NMF::predict(NMF::nmf(
              x.nmf.samp, rank = nk[k], method = "lee", seed = seed)),
            hcAEucl = hcAEucl(x.rest[ind.new, ], nk[k]),
            hcDianaEucl = hcDianaEucl(x.rest[ind.new, ], nk[k]),
            kmEucl = stats::kmeans(
              stats::dist(x.rest[ind.new, ]), nk[k])$cluster,
            kmSpear = stats::kmeans(
              spearman_dist(x.rest[ind.new, ]), nk[k])$cluster,
            pamEucl = cluster::pam(
              stats::dist(x.rest[ind.new, ]), nk[k], cluster.only = TRUE),
            pamSpear = cluster::pam(
              spearman_dist(x.rest[ind.new, ]), nk[k], cluster.only = TRUE),
            apEucl = stats::setNames(dplyr::dense_rank(suppressWarnings(
              apcluster::apclusterK(apcluster::negDistMat, x.rest[ind.new, ],
                                    nk[k], verbose = FALSE)@idx)),
              rownames(x.rest[ind.new, ])),
            scRbf = stats::setNames(kernlab::specc(x.rest[ind.new, ], nk[k],
                                                   kernel = "rbfdot")@.Data,
                                    rownames(x.rest[ind.new, ])),
            gmmBIC = mclust::Mclust(x.rest[ind.new, ], nk[k])$classification,
            biclust = blockcluster::cocluster(
              x.rest[ind.new, ], "continuous",
              nbcocluster = c(nk[k], nk[k]))@rowclass + 1
          )
        }
      }
      coclus[, , , k]
    } %>%
      unlist() %>% 
      array(dim = c(n, reps, nm, lnk),
            dimnames = list(samples, paste0("R", 1:reps), algorithms, nk))
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
#' prepare_data(x)
prepare_data <- function(data, min.sd = 1) {
  dat.out <- data %>%
    magrittr::extract(apply(., 1, function(x) !any(is.na(x))),
                      apply(., 2, function(x) stats::sd(x, na.rm = TRUE)) >
                        min.sd) %>%
    scale()
  return(dat.out)
}

#' Hierarchical clustering with Euclidean distance and Average linkage
#' @param d data matrix
#' @param k scalar indicating number of clusters to cut tree into
#' @noRd
hcAEucl <- function(d, k) {
  return(as.integer(stats::cutree(stats::hclust(
    stats::dist(d), method = "average"), k)))
}

#' Hierarchical clustering using DIvisive ANAlysis algorithm
#'
#' @inheritParams hcAEucl
#' @noRd
hcDianaEucl <- function(d, k) {
  return(as.integer(stats::cutree(cluster::diana(
    stats::dist(d), diss = TRUE), k)))
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
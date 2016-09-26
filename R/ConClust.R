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
#' \item{"hcSEucl": }{Hierarchical Clustering using Single linkage and Euclidean distance}
#' \item{"hcDianaEucl": }{Divisive Hierarchical Clustering using Euclidean distance}
#' \item{"kmEucl": }{K-Means using Euclidean distance}
#' \item{"kmSpear": }{K-Means using Spearman distance}
#' \item{"kmMI": }{K-Means using Mutual Information distance}
#' \item{"pamEucl": }{Partition Around Mediods using Euclidean distance}
#' \item{"pamSpear": }{Partition Around Mediods using Spearman distance}
#' \item{"pamMI": }{Partition Around Mediods using Mutual Information distance}
#' \item{"apEucl": }{Affinity Propagation using Euclidean distance}
#' \item{"scRbf": }{Spectral Clustering using Radial-Basis kernel function}
#' \item{"gmmBIC": }{Gaussian Mixture Model using Bayesian Information Criterion on EM algorithm}
#' \item{"biclust": }{Biclustering using a latent block model}
#' }
#'
#' The \code{min.sd} argument is used to filter the feature space for only highly variable
#' features. Only features with a standard deviation across all samples greater than
#' \code{min.sd} will be used.
#'
#' @param x data matrix; genes are rows and samples are columns
#' @param k number of clusters requested
#' @param pItem proportion of items to be used in subsampling within an algorithm
#' @param reps number of subsamples
#' @param method vector of clustering algorithms for performing consensus clustering. Must be
#' any number of the following: "nmfDiv", "nmfEucl", "hcAEucl", "hcDianaEucl", "kmEucl",
#' "kmSpear", "kmMI", "pamEucl", "pamSpear", "pamMI", "apEucl", "scRbf", "gmmBIC", "biclust".
#' See details.
#' @param seed random seed to use for NMF-based algorithms
#' @param seed.method seed to use to ensure each method operates on the same set of subsamples
#' @param min.sd minimum standard deviation threshold. See details.
#' @param dir directory where returned object will be saved at each iteration (as an RDS object).
#' No output file is saved if \code{file} is \code{NULL}.
#' @param fileName file name of the written object
#' @param time.saved logical; if \code{TRUE} (default), the date saved is appended
#' to the file name. Only applicable when \code{dir} is not \code{NULL}.
#' @return An array of dimension \code{nrow(x)} by \code{reps} by \code{length(methods)}
#' Each slice of the array is a matrix showing consensus clustering results for
#' algorithms (method). The matrices have a row for each sample, and a column for each
#' subsample. Each entry represents a class membership.
#' @author Derek Chiu, Aline Talhouk
#' @importFrom magrittr extract
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats dist hclust kmeans setNames
#' @import apcluster kernlab mclust blockcluster
#' @export
ConClust <- function(x, k, pItem = 0.8, reps = 1000, method = NULL,
                     seed = 123456, seed.method = 1, min.sd = 1, dir = NULL,
                     fileName = "ConClustOutput", time.saved = TRUE) {
  . <- NULL
  x.rest <- dataPrep(x, min.sd = min.sd)
  x.nmf <- x.rest %>%
    rbind(-.) %>%
    apply(2, function(x) ifelse(x < 0, 0, x))

  samples <- colnames(x.rest)
  n <- ncol(x.rest)
  n.new <- floor(n * pItem)
  if (is.null(method))
    method <- c("nmfDiv", "nmfEucl", "hcAEucl", "hcDianaEucl", "kmEucl",
                "kmSpear", "pamEucl", "pamSpear", "apEucl",
                "scRbf", "gmmBIC", "biclust")
  nm <- length(method)
  coclus <- array(NA, c(n, reps, nm),
                  dimnames = list(samples, paste0("R", 1:reps), method))
  pb <- txtProgressBar(min = 0, max = reps, style = 3)

  for (j in 1:nm) {
    set.seed(seed.method)
    for (i in 1:reps) {
      setTxtProgressBar(pb, i)
      ind.new <- sample(n, n.new, replace = F)
      if (any(c("nmfDiv", "nmfEucl") %in% method))
        x.nmf.samp <- x.nmf[!(apply(x.nmf[, ind.new], 1,
                                    function(x) all(x == 0))), ind.new]

      coclus[ind.new, i, j] <- switch(
        method[j],
        nmfDiv = NMF::predict(NMF::nmf(
          x.nmf.samp, rank = k, method = "brunet", seed = seed)),
        nmfEucl = NMF::predict(NMF::nmf(
          x.nmf.samp, rank = k, method = "lee", seed = seed)),
        hcAEucl = cutree(hclust(dist(t(x.rest[, ind.new])),
                                method = "average"), k),
        hcSEucl = cutree(hclust(dist(t(x.rest[, ind.new])),
                                method = "single"), k),
        hcDianaEucl = cutree(diana(euc(t(x.rest[, ind.new])),
                                   diss = TRUE), k),
        kmEucl = kmeans(euc(t(x.rest[, ind.new])),
                        k)$cluster,
        kmSpear = kmeans(spearman.dist(t(x.rest[, ind.new])),
                         k)$cluster,
        pamEucl = pam(euc(t(x.rest[, ind.new])), k,
                      cluster.only = TRUE),
        pamSpear = pam(spearman.dist(t(x.rest[, ind.new])), k,
                       cluster.only = TRUE),
        apEucl = setNames(dense_rank(suppressWarnings(
          apclusterK(negDistMat, t(x.rest[, ind.new]), k,
                     verbose = FALSE)@idx)),
          rownames(t(x.rest[, ind.new]))),
        scRbf = setNames(specc(t(x.rest[, ind.new]), k,
                               kernel = "rbfdot")@.Data,
                         rownames(t(x.rest[, ind.new]))),
        gmmBIC = Mclust(t(x.rest[, ind.new]), k)$classification,
        biclust = cocluster(as.matrix(x.rest[, ind.new]), "continuous",
                            nbcocluster = c(k, k))@colclass + 1
        )
      if(i %% 10==0) {
          saveRDS(coclus, paste0(dir, fileName,".rds"))
        }
    }
  }
  saveRDS(coclus, paste0(dir, fileName,".rds"))
  return(coclus)
}

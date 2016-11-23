#' Testing feature/custom-algs
#' @export
#' @examples 
#' data(hgsc)
#' dat <- t(hgsc[, -1])
#' 
#' # Custom distance function
#' manh = function(x) {
#'   stats::dist(x, method = "manhattan")
#' }
#' 
#' # Custom clustering algorithm
#' agnes <- function(d, k) {
#'   return(as.integer(stats::cutree(cluster::agnes(d, diss = TRUE), k)))
#' }
#' 
#' cc <- consensus_cluster2(dat, reps = 5, algorithms = c("hc", "agnes"),
#' distance = c("euclidean", "manh"))
#' str(cc)
consensus_cluster2 <- function(data, nk = 2:4, pItem = 0.8, reps = 1000,
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
  return(all.arr)
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
#' Testing feature/custom-algs
#' @export
#' @examples 
#' data(hgsc)
#' dat <- t(hgsc[, -1])
#' test <- consensus_cluster(dat, reps = 5)
#' test2 <- consensus_cluster2(dat, reps = 5)
#' str(test)
#' str(test2)
consensus_cluster2 <- function(data, nk = 2:4, pItem = 0.8, reps = 1000,
                               algorithms = NULL, nmf.method = c("brunet", "lee"),
                               distance = "euclidean",
                               parallel = NULL, ncores = NULL,
                               progress = TRUE, seed = 123456, seed.alg = 1,
                               min.sd = 1, save = FALSE, file.name = "CCOutput",
                               time.saved = FALSE) {
  x.rest <- prepare_data(data, min.sd = min.sd)
  if (is.null(algorithms))
    algorithms <- c("nmf", "hc", "diana", "km", "pam",
                    "ap", "sc", "gmm", "block")
  if ("nmf" %in% algorithms) {
    nmf.arr <- cluster_nmf(x.rest, nk, pItem, reps, nmf.method, seed, seed.alg)
  } else {
    nmf.arr <- NULL
  }
  
  oalgs <- c("ap", "sc", "gmm", "block")
  if (any(oalgs %in% algorithms)) {
    other.arr <- cluster_other(x.rest, nk, pItem, reps, oalgs, seed, seed.alg)
  } else {
    other.arr <- NULL
  }
  
  dalgs <- c("hc", "diana", "km", "pam")
  if (any(dalgs %in% algorithms)) {
    dist.arr <- cluster_dist(x.rest, nk, pItem, reps, dalgs, distance, seed, seed.alg)
  } else {
    dist.arr <- NULL
  }

  all.arr <- abind::abind(nmf.arr, dist.arr, other.arr, along = 3)
  return(all.arr)
}

#' Cluster NMF-based algorithms
#' @noRd
cluster_nmf <- function(data, nk, pItem, reps, nmf.method, seed, seed.alg) {
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
      }
    }
  }
  return(nmf.arr)
}

#' Cluster algorithms with dissimilarity specification
#' @noRd
cluster_dist <- function(data, nk, pItem, reps, dalgs, distance, seed, seed.alg) {
  n <- nrow(data)
  n.new <- floor(n * pItem)
  ld <- length(distance)
  ldist <- prod(length(dalgs), ld)
  lnk <- length(nk)
  dist.arr <- array(NA, c(n, reps, ldist, lnk),
                    dimnames = list(rownames(data),
                                    paste0("R", 1:reps),
                                    apply(expand.grid(Hmisc::capitalize(distance),
                                                      toupper(dalgs)),
                                          1, function(x) paste0(x[2], "_", x[1])),
                                    nk))
  for (k in 1:lnk) {
    set.seed(seed.alg)
    for (i in 1:reps) {
      ind.new <- sample(n, n.new, replace = FALSE)
      dists <- distances(data[ind.new, ], distance)
      for (j in 1:length(dalgs)) {
        for (d in 1:ld) {
          dist.arr[ind.new, i, (j - 1) * ld + d, k] <- get(dalgs[j])(dists[[d]], 4)
        }
      }
    }
  }
  return(dist.arr)
}

#' Cluster other algorithms
#' @noRd
cluster_other <- function(data, nk, pItem, reps, oalgs, seed, seed.alg) {
  n <- nrow(data)
  n.new <- floor(n * pItem)
  lalgs <- length(oalgs)
  lnk <- length(nk)
  other.arr <- array(NA, c(n, reps, lalgs, lnk),
                     dimnames = list(rownames(data),
                                     paste0("R", 1:reps),
                                     toupper(oalgs),
                                     nk))
  for (k in 1:lnk) {
    for (j in 1:length(oalgs)) {
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
      }
    }
  }
  return(other.arr)
}

#' @noRd
myDistFunc = function(x) {
  stats::dist(x, method = "manhattan")
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

# algorithms <- c("nmf", "hc", "diana", "km", "pam", "ap", "sc", "gmm", "block")
# nalgs <- "nmf"
# dalgs <- c("hc", "diana", "km", "pam")
# oalgs <- c("ap", "sc", "gmm", "block")
# distance <- c("eucl", "spear", "myDistFunc")
# linkage <- "average"
# nmf.method <- c("brunet", "lee")
# nnmf <- length(nmf.method)
# methods <- c("NMF_Div", "NMF_Eucl", "HC_Eucl", "HC_Spear",
#              "DIANA_Eucl", "DIANA_Spear", "KM_Eucl", "KM_Spear",
#              "PAM_Eucl", "PAM_Spear", "AP", "SC", "GMM", "BLOCK")
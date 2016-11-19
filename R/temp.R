
# Setup -------------------------------------------------------------------
myDistFunc = function(x) {
  stats::dist(x, method = "manhattan")
}

hc <- function(d, k, method = "average") {
  return(as.integer(stats::cutree(stats::hclust(d, method = method), k)))
}

diana <- function(d, k) {
  return(as.integer(stats::cutree(cluster::diana(d, diss = TRUE), k)))
}

km <- function(d, k) {
  return(as.integer(stats::kmeans(d, k)$cluster))
}

pam <- function(d, k) {
  return(as.integer(cluster::pam(d, k, cluster.only = TRUE)))
}

algorithms <- c("nmf", "hc", "diana", "km", "pam", "ap", "sc", "gmm", "block")
nalgs <- "nmf"
dalgs <- c("hc", "diana", "km", "pam")
oalgs <- c("ap", "sc", "gmm", "block")
distance <- c("eucl", "spear", "myDistFunc")
linkage <- "average"
nmf.method <- c("brunet", "lee")
nnmf <- length(nmf.method)
methods <- c("NMF_Div", "NMF_Eucl", "HC_Eucl", "HC_Spear",
             "DIANA_Eucl", "DIANA_Spear", "KM_Eucl", "KM_Spear",
             "PAM_Eucl", "PAM_Spear", "AP", "SC", "GMM", "BLOCK")


# NMF ---------------------------------------------------------------------
#' Cluster NMF-based algorithms
cluster_nmf <- function(data, reps, nmf.method, nk, seed, seed.alg) {
  x.nmf <- data %>%
    cbind(-.) %>%
    apply(2, function(x) ifelse(x < 0, 0, x))
  n <- nrow(data)
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




# Distance ----------------------------------------------------------------
#' Cluster algorithms with dissimilarity specification
cluster_dist <- function(data, reps, dalgs, distance, nk, seed, seed.alg) {
  n <- nrow(data)
  ld <- length(distance)
  ldist <- prod(length(dalgs), ld)
  lnk <- length(nk)
  dist.arr <- array(NA, c(n, reps, ldist, lnk),
                    dimnames = list(samples,
                                    paste0("R", 1:reps),
                                    apply(expand.grid(Hmisc::capitalize(distance),
                                                      toupper(dalgs)),
                                          1, function(x) paste0(x[2], "_", x[1])),
                                    nk))
  for (k in 1:lnk) {
    set.seed(seed.alg)
    for (i in 1:reps) {
      ind.new <- sample(n, n.new, replace = FALSE)
      dists <- distances(x.rest[ind.new, ], distance)
      for (j in 1:length(dalgs)) {
        for (d in 1:ld) {
          dist.arr[ind.new, i, (j - 1) * ld + d, k] <- get(dalgs[j])(dists[[d]], 4)
        }
      }
    }
  }
  return(dist.arr)
}


# Other -------------------------------------------------------------------
#' Cluster other algorithms

cluster_other <- function(data, reps, oalgs, nk, seed, seed.alg) {
  n <- nrow(data)
  lalgs <- length(oalgs)
  lnk <- length(nk)
  other.arr <- array(NA, c(n, reps, lalgs, lnk),
                     dimnames = list(samples,
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
                   apcluster::apclusterK(apcluster::negDistMat, x.rest[ind.new, ],
                                         nk[k], verbose = FALSE)@idx)),
                   rownames(x.rest[ind.new, ])),
                 sc = stats::setNames(kernlab::specc(x.rest[ind.new, ], nk[k],
                                                     kernel = "rbfdot")@.Data,
                                      rownames(x.rest[ind.new, ])),
                 gmm = mclust::Mclust(x.rest[ind.new, ], nk[k])$classification,
                 block = blockcluster::cocluster(
                   x.rest[ind.new, ], "continuous",
                   nbcocluster = c(nk[k], nk[k]))@rowclass + 1)
      }
    }
  }
  return(other.arr)
}


# Example -----------------------------------------------------------------
nmf.arr <- cluster_nmf(x.rest, reps, nmf.method, nk, seed, seed.alg)
dist.arr <- cluster_dist(x.rest, reps, dalgs, distance, nk, seed, seed.alg)
other.arr <- cluster_other(x.rest, reps, oalgs, nk, seed, seed.alg)
all.arr <- abind::abind(nmf.arr, dist.arr, other.arr, along = 3)

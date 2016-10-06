#' Title function for computing compactness validity index for a clustering result
#'
#' @param X : a dataset, rows are observations, columns are variables
#' @param labels : cluster labels from a clustering result(a vector)
#'
#' @return compactness score
#' @export
#'
#' @examples
#' data("FGD")
#' data("E")
#' labels<-E[,1]
#' valid_compactness(FGD,labels)
valid_compactness <- function(X, labels) {
  assertthat::assert_that(is.data.frame(X), length(labels) == nrow(X))
  n <- length(labels)
  C <- sort(unique(labels))
  k <- length(C)
  cp <- 0
  for (i in 1:k) {
    ind <- which(labels == C[i])
    nk <- length(ind)
    if (nk <= 1) {
      cp <- cp + 0
    } else{
      sum_d <- 0
      sum_d <- sum(dist(X[ind, ], method = "euclidean"))
      cp <- cp + (nk * (sum_d / (nk * (nk - 1) / 2)))
    }
  }
  return(cp / n)
}

#' Title function for computing Davies-Bouldin index and Dunn index
#'
#' @param X : a data set whose rows are observations, columns are variables
#' @param labels : cluster labels from a clustering result(a vector)
#'
#' @return DB: Davies-Bouldin score, Dunn: Dunn score
#' @export
#'
#' @examples
#' data("FGT")
#' data("FGD")
#' valid_DbDunn(FGD,FGT)
valid_DbDunn <- function(X, labels) {
  if (is.data.frame(labels)) {
    assertthat::assert_that(nrow(X) == nrow(labels))
  } else if (is.vector(labels)) {
    assertthat::assert_that(nrow(X) == length(labels))
  } else{
    stop(
      "labels is neither a vector nor a data frame with length equal to number of observations in the data set."
    )
  }
  assertthat::assert_that(is.data.frame(X))
  nrow <- nrow(X)
  nc <- ncol(X)
  k <- max(labels)
  temp <- valid_sumsqures(X, labels, k)
  st <- temp$Tot
  sw <- temp$W
  sb <- temp$B
  cintra <- temp$Sintra
  cinter <- temp$Sinter
  R <- zeros(k)
  dbs <- zeros(1, k)
  for (i in 1:k) {
    for (j in (i + 1):k) {
      if (j <= k) {
        if ((cinter[i, j] == 0)) {
          R[i, j] <- 0
        } else{
          R[i, j] <- (cintra[i] + cintra[j]) / cinter[i, j]
        }
      }
    }
    dbs[1, i] <- max(R[i, ])
  }
  DB <- mean(dbs[1, 1:k - 1])
  dbs <- max(cintra)
  R <- cinter / dbs
  for (i in 1:(k - 1)) {
    S <- R[i, (i + 1):k]
    dbs[i] <- min(S)
  }
  return(list(DB = DB, Dunn = min(dbs)))
}

#' Title function to compute within group, between group, and total sum of squares and cross-products
#'
#' @param data : a matrix with each column representing a variable
#' @param labels : a vector indicating class labels
#' @param k : number of clusters
#'
#' @return W: within group sum of squares and cross-products; B: between group sum of squares and cross-products;T: total sum of squares and cross-products;
#'          Sintra & Sinter: centroid diameter and linkage distance
#' @export
#'
#' @examples
#' data("FGT")
#' data("FGD")
#' valid_sumsqures(data=FGD,labels=FGT,k=4)
valid_sumsqures <- function(data, labels, k) {
  assertthat::assert_that(is.vector(labels) || is.data.frame(labels))
  if (is.data.frame(labels)) {
    assertthat::assert_that(nrow(data) == nrow(labels))
  } else{
    assertthat::assert_that(nrow(data) == length(labels))
  }
  assertthat::assert_that(is_pos_int(k))
  if (is.data.frame(data)) {
    data <- as.matrix(data)
  }
  ncase <- nrow(data)
  m <- ncol(data)
  Dm <- t(colMeans(data))
  Dm <- data - Dm[ones(ncase, 1), ]
  Tot <- t(Dm) %*% Dm
  W <- matrix(rep(0, nrow(Tot) * ncol(Tot)), nrow = nrow(Tot))
  Dm <- matrix(rep(0, k * m), nrow = k)
  Sintra <- matrix(rep(0, k), nrow = 1)
  for (i in 1:k) {
    if (k > 1) {
      Cindex <- which(labels == i)
    } else{
      Cindex <- 1:ncase
    }
    nk <- length(Cindex)
    if (nk > 1) {
      dataC <- data[Cindex, ]
      m <- colMeans(dataC)
      Dm[i, ] <- m
      dataC <- dataC - repmat(m, nk, 1)
      W <- W + t(dataC) %*% dataC
      dataC <- rowSums(dataC ^ 2)
      Sintra[i] <- mean(sqrt(dataC))
    }
  }
  B <- Tot - W
  Sinter <- matrix(rep(0, k ^ 2), nrow = k)
  if (k > 1) {
    for (i in 1:k) {
      for (j in (i + 1):k) {
        if (j <= k) {
          m <- abs(Dm[i, ] - Dm[j, ])
          Sinter[i, j] <- sqrt(sum(m ^ 2))
          Sinter[j, i] <- Sinter[i, j]
        }
      }
    }
  }
  return(list(
    Tot = Tot,
    W = W,
    B = B,
    Sintra = Sintra,
    Sinter = Sinter
  ))
}
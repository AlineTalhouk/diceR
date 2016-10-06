#' Internal validity indices
#' 
#' \strong{I}nternal \strong{v}alidity indices are agnostic to any reference
#' clustering results. They calculate cluster performance on the basis of 
#' compactness and separability.
#' 
#' \code{iv_compactness} computes the compactness validity index for a
#' clustering result
#'
#' @param x a dataset with rows as observations, columns sariables
#' @param labels a vector of cluster labels from a clustering result
#'
#' @return \code{iv_compactness} returns the compactness score
#' @name internal_validity
#' @export
#'
#' @examples
#' data("FGD", "E", "FGT", "FGD")
#' iv_compactness(FGD, E[, 1])
#' iv_db_dunn(FGD, FGT)
#' iv_sumsq(data = FGD, labels = FGT, k = 4)
iv_compactness <- function(x, labels) {
  assertthat::assert_that(is.data.frame(x), length(labels) == nrow(x))
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
      sum_d <- sum(dist(x[ind, ], method = "euclidean"))
      cp <- cp + (nk * (sum_d / (nk * (nk - 1) / 2)))
    }
  }
  return(cp / n)
}

#' @details \code{iv_db_dunn} computes the Davies-Bouldin index and Dunn index
#' @return \code{iv_db_dunn} returns a list with elements
#'   \item{DB}{Davies-Bouldin score}
#'   \item{Dunn}{Dunn score}
#' @rdname internal_validity
#' @export
iv_db_dunn <- function(x, labels) {
  if (is.data.frame(labels)) {
    assertthat::assert_that(nrow(x) == nrow(labels))
  } else if (is.vector(labels)) {
    assertthat::assert_that(nrow(x) == length(labels))
  } else{
    stop(
      "labels is neither a vector nor a data frame with length equal to number of observations in the data set."
    )
  }
  assertthat::assert_that(is.data.frame(x))
  nrow <- nrow(x)
  nc <- ncol(x)
  k <- max(labels)
  temp <- iv_sumsq(x, labels, k)
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

#' @details \code{iv_sumsq} computes the within group, between group, and total
#'   sum of squares and cross-products
#' @param k number of clusters
#'
#' @return \code{iv_sumsq} returns a list with elements 
#'   \item{W}{within group sum of squares and cross-products}
#'   \item{B}{between group sum of squares and cross-products}
#'   \item{T}{total sum of squares and cross-products}
#'   \item{Sintra}{centroid diameter}
#'   \item{Sinter}{linkage distance}          
#' @rdname internal_validity
#' @export
iv_sumsq <- function(x, labels, k) {
  assertthat::assert_that(is.vector(labels) || is.data.frame(labels))
  if (is.data.frame(labels)) {
    assertthat::assert_that(nrow(x) == nrow(labels))
  } else{
    assertthat::assert_that(nrow(x) == length(labels))
  }
  assertthat::assert_that(is_pos_int(k))
  if (is.data.frame(x)) {
    data <- as.matrix(x)
  }
  ncase <- nrow(x)
  m <- ncol(x)
  Dm <- t(colMeans(data))
  Dm <- x - Dm[ones(ncase, 1), ]
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
      dataC <- x[Cindex, ]
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
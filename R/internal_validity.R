#' Internal validity indices
#' 
#' \strong{I}nternal \strong{v}alidity indices are agnostic to any reference 
#' clustering results. They calculate cluster performance on the basis of 
#' compactness and separability.
#' 
#' \code{iv_compactness} computes the compactness validity index for a 
#' clustering result
#' 
#' @param data a dataset with rows as observations, columns as variables
#' @param labels a vector of cluster labels from a clustering result
#' @references MATLAB functions \code{valid_compactness}, \code{valid_DbDunn}, 
#'   \code{valid_sumsqures} by Simon Garrett in LinkCluE
#' @references Calinski, T., Harabasz, J. (1974), A dendrite method for cluster
#'   analysis, "Communications in Statistics", vol. 3, 1-27.
#' @return \code{iv_compactness} returns the compactness score
#' @name internal_validity
#' @export
#' 
#' @examples
#' set.seed(1)
#' E <- matrix(rep(sample(1:4, 1000, replace = TRUE)), nrow = 100, byrow =
#'               FALSE)
#' set.seed(1)
#' dat <- as.data.frame(matrix(runif(1000, -10, 10), nrow = 100, byrow = FALSE))
#' iv_compactness(dat, E[, 1])
#' iv_db_dunn(dat, E[, 1])
#' iv_chi(dat, E[, 1], centrotypes = "centroids")
iv_compactness <- function(data, labels) {
  assertthat::assert_that(is.data.frame(data), length(labels) == nrow(data))
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
      sum_d <- sum(dist(data[ind, ], method = "euclidean"))
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
iv_db_dunn <- function(data, labels) {
  if (is.data.frame(labels)) {
    assertthat::assert_that(nrow(data) == nrow(labels))
  } else if (is.vector(labels)) {
    assertthat::assert_that(nrow(data) == length(labels))
  } else{
    stop(
      "labels is neither a vector nor a data frame with length equal to number of observations in the data set."
    )
  }
  assertthat::assert_that(is.data.frame(data))
  nrow <- nrow(data)
  nc <- ncol(data)
  k <- max(labels)
  temp <- sumsq(data, labels, k)
  st <- temp$Tot
  sw <- temp$W
  sb <- temp$B
  cintra <- temp$Sintra
  cinter <- temp$Sinter
  R <- pracma::zeros(k)
  dbs <- pracma::zeros(1, k)
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

#' Sum of squares
#' @param k number of clusters
#' @return \code{sumsq} returns a list with elements 
#'   \item{W}{within group sum of squares and cross-products}
#'   \item{B}{between group sum of squares and cross-products}
#'   \item{T}{total sum of squares and cross-products}
#'   \item{Sintra}{centroid diameter}
#'   \item{Sinter}{linkage distance}   
#' @noRd       
sumsq <- function(data, labels, k) {
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
  Dm <- data - Dm[pracma::ones(ncase, 1), ]
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
      dataC <- dataC - pracma::repmat(m, nk, 1)
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

#' @details \code{iv_chi} computes the Calinski-Harabasz Index, or pseudo
#'   F-statistic. The function is a wrapper for
#'   \code{\link[clusterSim]{index.G1}}
#' @param ... additional arguments to \code{\link[clusterSim]{index.G1}}
#' @return \code{iv_chi} returns the CHI
#' @rdname internal_validity
#' @export
iv_chi <- function(data, ...) {
  clusterSim::index.G1(data, ...)
}

#' Proportion of Ambiguous Clustering
#'
#' Given a consensus matrix, returns the proportion of ambiguous clusters (PAC).
#' This is a robust way to assess clustering performance.
#'
#' Since a consensus matrix is symmetric, we only look at its lower (or upper)
#' triangular matrix. The proportion of entries strictly between \code{lower}
#' and \code{upper} is the PAC. In a perfect clustering, the consensus matrix
#' would consist of only 0s and 1s, and the PAC assessed on the (0, 1) interval
#' would have a perfect score of 0. Using a (0.1, 0.9) interval for defining
#' ambiguity is common as well.
#' 
#' The PAC is not, strictly speaking, an internal validity index. Originally
#' used to choose the optimal number of clusters, here we use it to assess
#' cluster stability. However, PAC is still agnostic any gold standard
#' clustering result so we use it like an internal validity index.
#'
#' @param cm consensus matrix. Should be symmetric and values between 0
#'   and 1.
#' @param lower the lower bound that determines what is ambiguous
#' @param upper the upper bound that determines what is ambiguous
#' @return the PAC is a score used in clustering performance. The lower it is
#'   the better, because we want minimal ambiguity amongst the consensus.
#' @author Derek Chiu
#' @references Senbabaoglu, Y., Michailidis, G., & Li, J. Z. (2014). Critical
#'   limitations of consensus clustering in class discovery. Scientific reports,
#'   4.
#' @export
#' @examples
#' set.seed(1)
#' x <- replicate(100, rbinom(100, 4, 0.2))
#' y <- consensus_matrix(x)
#' PAC(y, lower = 0.05, upper = 0.95)
PAC <- function(cm, lower = 0, upper = 1) {
  pac <- cm %>%
    extract(lower.tri(.)) %>%
    extract(and(is_greater_than(., lower), is_less_than(., upper))) %>%
    length() %>%
    divide_by(., length(cm[lower.tri(cm)]))
  return(pac)
}
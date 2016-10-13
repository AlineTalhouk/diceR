#' Generate clusters using distance criterion
#' 
#' The function generates clusters using the distance criterion from a
#' hierarchical cluster tree, equivalent to \code{cluster} in MATLAB
#'
#' @param Z matrix representing a hierarchical cluster tree
#' @param maxclust desired number of clusters
#' @references MATLAB function cluster
#' @return clustering result
#' @export
#'
#' @examples
#' data("E_LCE")
#' tree <- linkage(stod(asrs(E_LCE, 0.8)), "complete")
#' cluster_tree <- dist_cluster(tree, 4)
dist_cluster <- function(Z, maxclust) {
  assertthat::assert_that(is.matrix(Z))
  assertthat::assert_that(is.numeric(Z))
  assertthat::assert_that(is_pos_int(maxclust))
  m <- nrow(Z) + 1
  resultT <- pracma::zeros(m, length(maxclust))
  for (j in 1:length(maxclust)) {
    if (m <= maxclust[j]) {
      resultT[, j] <- 1:m
    } else if (maxclust[j] == 1) {
      resultT[, j] <- rep(1, m)
    } else{
      clsnum <- 1
      for (k in (m + 1 - maxclust[j]):(m - 1)) {
        i <- Z[k, 1]
        if (i <= m) {
          resultT[i, j] <- clsnum
          clsnum <- clsnum + 1
        } else if (i < (2 * m - maxclust[j] + 1)) {
          resultT[, j] <- clusternum(Z, resultT[, j], i - m, clsnum)
          clsnum <- clsnum + 1
        }
        i <- Z[k, 2]
        if (i <= m) {
          resultT[i, j] <- clsnum
          clsnum <- clsnum + 1
        } else if (i < (2 * m - maxclust[j] + 1)) {
          resultT[, j] <- clusternum(Z, resultT[, j], i - m, clsnum)
          clsnum <- clsnum + 1
        }
      }
    }
  }
  return(resultT)
}

#' Assign leaves under cluster c to c, equivalent to MATLAB function clusternum
#'
#' @param X original hierarchical cluster tree
#' @param resultT cluster tree generated in cluster
#' @param k level in the tree
#' @param c cluster to be assigned to the leaf
#' @noRd
clusternum <- function(X, resultT, k, c) {
  assertthat::assert_that(is.numeric(X))
  assertthat::assert_that(is.numeric(resultT))
  m <- nrow(X) + 1
  children <- NULL
  while (length(k) != 0) {
    children <- X[k, 1:2]
    resultT[children[which(children <= m)]] <- c
    k <- children[which(children > m)] - m
  }
  return(resultT)
}
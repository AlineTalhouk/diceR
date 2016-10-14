#' Relabel clusters
#' 
#' Relabel clusters in the ensemble \code{E}
#'
#' @param E N by M matrix of cluster ensemble
#'
#' @return A list with elements
#' \item{newE}{N by M matrix of relabelled ensemble}
#' \item{no_allcl}{total number of clusters in the ensemble}
#' @references MATLAB function relabelCl by Simon Garrett in LinkCluE package
#' @author Johnson Liu
#' @export
#'
#' @examples
#' set.seed(1)
#' E<-matrix(rep(sample(1:4,1000,replace = TRUE)),nrow=100,byrow=FALSE)
#' E_relabelled <- relabel_clusters(E)
relabel_clusters <- function(E) {
  assertthat::assert_that(is.matrix(E))
  assertthat::assert_that(is.numeric(E))
  N <- nrow(E)
  M <- ncol(E)
  for (i in 1:N) {
    for (j in 1:M) {
      if (!is_pos_int(E[i, j])) {
        stop("Error: one of the entries in the input matrix is not a positive integer.")
      }
    }
  }
  newE <- matrix(rep(0, N * M), nrow = N)
  ucl <- sort(unique(E[, 1]))
  if (max(E[, 1] != length(ucl)) == 1) {
    for (j in 1:length(ucl)) {
      newE[which(E[, 1] == ucl[j]), 1] <- j
    }
  }
  for (i in 2:M) {
    ucl <- sort(unique(E[, i]))
    prevCl <- length(sort(unique(c(newE[, 1:(i - 1)]))))
    for (j in 1:length(ucl)) {
      newE[E[, i] == ucl[j], i] <- prevCl + j
    }
  }
  no_allcl <- max(max(newE))
  return(list(no_allcl = no_allcl, newE = newE))
}
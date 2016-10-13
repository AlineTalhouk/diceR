#' Weigh clusters
#' 
#' Compute weight for each pair of clusters using their shared members (Jaccard
#' coefficient)
#' 
#' @param E N by M matrix of cluster ensembles
#' @references MATALAB function weightCl by Simon Garrett in package LinkCluE   
#' @return a p by p weighted cluster matrix where p denotes number of classes
#' @author Johnson Liu
#' @export
#' 
#' @examples
#' data("E_LCE")
#' wl <- weigh_clusters(E_LCE)
weigh_clusters <- function(E) {
  assertthat::assert_that(is.matrix(E), is.numeric(E))
  for (i in 1:nrow(E)) {
    for (j in 1:ncol(E)) {
      if (!is_pos_int(E[i, j])) {
        stop("Error in relabel_clusters: one of the entries in the input matrix is not a positive integer.")
      }
    }
  }
  N = nrow(E)
  no_allcl <- max(max(E))
  pc <- matrix(rep(0, N * no_allcl), nrow = N)
  for (i in 1:N) {
    pc[i, E[i, ]] <- 1
  }
  
  wcl <- matrix(rep(0, no_allcl ^ 2), nrow = no_allcl)
  for (i in 1:(no_allcl - 1)) {
    for (j in (i + 1):no_allcl) {
      tmp <- sum(as.numeric((pc[, i] + pc[, j])) > 0)
      if (tmp > 0) {
        wcl[i, j] <- sum(as.numeric((pc[, i] + pc[, j])) == 2) / tmp
      }
    }
  }
  wcl <- wcl + t(wcl)
  return(wcl)
}
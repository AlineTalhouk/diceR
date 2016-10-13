#' Title function to take many clustering results and create a cluster ensemble
#'
#' @param E : matrix of clustering resutls 
#' @param dcCTS : decay constant for CTS matrix
#' @param dcSRS : decay constant for SRS matrix
#' @param dcASRS : decay constant for ASRS matrix
#' @param R : number of repetitions for SRS matrix
#' @param is.relabelled : boolean representing whether input E is relabelled
#'
#' @return a list containing final clustering result using k modes, majority vote, and hierarchical clustering(average, single, and complete linkages)
#' @export
#'
#' @examples
link_clust <-
  function(E,
           dcCTS = 0.8,
           dcSRS = 0.8,
           dcASRS = 0.8,
           R = 10,
           is.relabelled = TRUE) {
    assertthat::assert_that(is.array(E) || is.matrix(E))
    assertthat::assert_that(dcCTS >= 0 && dcCTS <= 1)
    assertthat::assert_that(dcASRS >= 0 && dcASRS <= 1)
    assertthat::assert_that(dcSRS >= 0 && dcSRS <= 1)
    assertthat::assert_that(is_pos_int(R))
    cluster_maj_vote <- NULL
    cluster_k_modes <- NULL
    cluster_hc <- NULL
    mNames <- c("cts", "srs", "asrs")
    dcs <- c(dcCTS, dcSRS, dcASRS)
    if (!is.na(dim(E)[3])) {
      maj_vote <- majority_voting(E = E, is.relabelled = is.relabelled)
      cluster_k_modes <- k_modes(E = E, is.relabelled = is.relabelled)
      flat_E <- E
      dim(flat_E) <- c(dim(E)[1], dim(E)[2] * dim(E)[3])
      if (is.relabelled == FALSE) {
        E_f <- apply(flat_E[,-1], 2, function(x) {
          relabel_class(x, flat_E[, 1])
        })
        flat_E <- cbind(flat_E[, 1], E_f)
      }
      if (anyNA(flat_E)) {
        flat_E <- t(apply(flat_E, 1, function(x) {
          x[which(is.na(x))] <- names(which.max(table(x)))
          return(x)
        }))
      }
      flat_E<-apply(flat_E,c(1,2),as.numeric)
      for (i in 1:length(mNames)) {
        if (mNames[i] == "srs") {
          S <- do.call(mNames[i], list(flat_E, dcs[i], R))
        } else{
          S <- do.call(mNames[i], list(flat_E, dcs[i]))
        }
        cluster_hc <- cbind(cluster_hc, clHC(S, K = max(flat_E)))
      }
    }
    return(list(cluster_hc=cluster_hc,cluster_k_modes=cluster_k_modes,cluster_maj_vote=cluster_maj_vote))
  }


#' Title: function to perform final clustering using hierarchical algorithms (single linkage-SL, complete linkage-CL, average linkage-AL)
#'
#' @param S : N by N similarity matrix
#' @param K : preferred number of clusters
#'
#' @return CR: an N by 3 matrix of clustering results from SL, CL, and AL
#' @export
#'
#' @examples
clHC <- function(S, K) {
  CR <- NULL
  d <- diceR::stod(S)
  Z <- diceR::linkage(d, "single")
  CR <- cbind(CR, diceR::dist_cluster(Z, K))
  Z <- diceR::linkage(d, "complete")
  CR <- cbind(CR, diceR::dist_cluster(Z, K))
  Z <- diceR::linkage(d, "average")
  CR <- cbind(CR, diceR::dist_cluster(Z, K))
  return(CR)
}
#' Title wrapper function to generate CTS, SRS, and ASRS matrices. If input has NA, majority voting is used to eliminate NA.
#'
#' @param E : matrix of clustering resutls (may be 3D or 2D)
#' @param dcCTS : decay constant for CTS matrix
#' @param dcSRS : decay constant for SRS matrix
#' @param dcASRS : decay constant for ASRS matrix
#' @param R : number of repetitions for SRS matrix
#' @param is.relabelled : boolean representing whether input E is relabelled
#' @author Johnson Liu
#' @return a list containing the CTS, SRS, and ASRS matrix
#' @export
#'
#' @examples
#' set.seed(1)
#' E <- matrix(rep(sample(1:4, 1000, replace = TRUE)), nrow = 100, byrow = FALSE)
#' link_clust(E=E,dcCTS=0.8,dcSRS=0.8,dcASRS=0.8,R=10,is.relabelled=FALSE)
#' data(hgsc)
#' x <- ConClust(hgsc[1:50], k = 4, reps = 4, method = c("nmfEucl", "hcAEucl",
#' "hcDianaEucl"), save = FALSE)
#' y <- link_clust(E=x,dcCTS=0.8,dcSRS=0.8,dcASRS=0.8,R=10)
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
    mNames <- c("cts", "srs", "asrs")
    dcs <- c(dcCTS, dcSRS, dcASRS)
    if (!is.na(dim(E)[3])) {
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
      flat_E <- apply(flat_E, c(1, 2), as.numeric)
      return(list(
        cts = cts(E = flat_E, dc = dcCTS),
        srs = srs(E = flat_E, dc = dcSRS, R = R),
        asrs = asrs(E = flat_E, dc = dcASRS)
      ))
    } else{
      if (sum(!complete.cases(E)) == 0) {
        return(list(
          cts = cts(E = E, dc = dcCTS),
          srs = srs(E = E, dc = dcSRS, R = R),
          asrs = asrs(E = E, dc = dcASRS)
        ))
      } else{
        E[which(is.na(E))] <- names(which.max(table(E)))
        E <- apply(E, c(1, 2), as.numeric)
        return(list(
          cts = cts(E = E, dc = dcCTS),
          srs = srs(E = E, dc = dcSRS, R = R),
          asrs = asrs(E = E, dc = dcASRS)
        ))
      }
    }
  }



#' Title: function to perform final clustering using hierarchical algorithms (single linkage-SL, complete linkage-CL, average linkage-AL)
#'
#' @param S : N by N similarity matrix
#' @param K : preferred number of clusters
#'
#' @return CR: an N by 3 matrix of clustering results from SL, CL, and AL
#' @references MATLAB function clHC in LinkCluE by Simon Garrett
#' @export
#'
#' @examples
#' set.seed(1)
#' E<-matrix(rep(sample(1:4,1000,replace = TRUE)),nrow=100,byrow=FALSE)
#' clHC(S=cts(E,0.8),K=4)
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
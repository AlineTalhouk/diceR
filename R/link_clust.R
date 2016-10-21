#' Generate CTS, SRS, and ASRS matrices.
#' 
#' If input has \code{NA}, majority voting is used to eliminate the missing
#' values.
#' 
#' @param E matrix of clustering results (may be 3D or 2D)
#' @param dcCTS decay constant for CTS matrix
#' @param dcSRS decay constant for SRS matrix
#' @param dcASRS decay constant for ASRS matrix
#' @param R number of repetitions for SRS matrix
#' @param is.relabelled boolean representing whether input E is relabelled
#' @author Johnson Liu
#' @return a list containing the CTS, SRS, and ASRS matrix
#' @importFrom stats complete.cases
#' @export
#' 
#' @examples
#' set.seed(1)
#' E <- matrix(rep(sample(1:4, 1000, replace = TRUE)), nrow = 100, byrow = FALSE)
#' link_clust(E = E, dcCTS = 0.8, dcSRS = 0.8, dcASRS = 0.8, R = 10,
#'            is.relabelled = FALSE)
#' data(hgsc)
#' dat <- t(hgsc[, -1])
#' x <- ConClust(dat[1:100, 1:50], k = 4, reps = 4,
#'               method = c("nmfEucl", "hcAEucl", "hcDianaEucl"), save = FALSE)
#' y <- link_clust(E = x, dcCTS = 0.8, dcSRS = 0.8, dcASRS = 0.8, R = 10)
link_clust <- function(E, dcCTS = 0.8, dcSRS = 0.8, dcASRS = 0.8, R = 10,
                       is.relabelled = TRUE) {
  assertthat::assert_that(is.array(E) || is.matrix(E), dcCTS >= 0 && dcCTS <= 1,
                          dcASRS >= 0 && dcASRS <= 1, dcSRS >= 0 && dcSRS <= 1,
                          is_pos_int(R))
  dcs <- c(dcCTS, dcSRS, dcASRS)
  if (!is.na(dim(E)[3])) {
    flat_E <- E
    dim(flat_E) <- c(dim(E)[1], dim(E)[2] * dim(E)[3])
    if (is.relabelled == FALSE) {
      E_f <- apply(flat_E[, -1], 2, function(x) {
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
  } else {
    if (sum(!complete.cases(E)) == 0) {
      return(list(
        cts = cts(E = E, dc = dcCTS),
        srs = srs(E = E, dc = dcSRS, R = R),
        asrs = asrs(E = E, dc = dcASRS)
      ))
    } else {
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
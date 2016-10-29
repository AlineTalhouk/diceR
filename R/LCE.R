#' Linkage Clustering Ensemble
#' 
#' Generate CTS, SRS, and ASRS similarity matrices.
#' 
#' @param E is an array of clustering results. Missing entries are imputed using
#'   \code{imputeMissing}
#' @param data original data matrix with rows as samples, columns as variables
#' @param k requested number of clusters
#' @param dcCTS decay constant for CTS matrix
#' @param dcSRS decay constant for SRS matrix
#' @param dcASRS decay constant for ASRS matrix
#' @param R number of repetitions for SRS matrix
#' @param sim.mat vector of similarity matrices to compute. Choices are "cts",
#'   "srs", "asrs".
#' @author Johnson Liu
#' @return a vector containing the cluster assignment from either the CTS, SRS,
#'   or ASRS similarity matrices
#' @export
#' @examples
#' data(hgsc)
#' dat <- t(hgsc[, -1])[1:100, 1:50]
#' x <- ConClust(dat, nc = 4, reps = 4,
#'               method = c("nmfEucl", "hcAEucl", "hcDianaEucl"), save = FALSE)
#' LCE(E = x, data = dat, k = 4, dcCTS = 0.8, dcSRS = 0.8, dcASRS = 0.8, R = 10,
#' sim.mat = "asrs")
LCE <- function(E, data, k, dcCTS = 0.8, dcSRS = 0.8, dcASRS = 0.8, R = 10,
                sim.mat = c("cts", "srs", "asrs")) {
  assertthat::assert_that(is.array(E), dcCTS >= 0 && dcCTS <= 1,
                          dcASRS >= 0 && dcASRS <= 1, dcSRS >= 0 && dcSRS <= 1,
                          is_pos_int(R))
  E2 <- imputeMissing(E, data, imputeALL = TRUE)$E_imputed2[, , 1]
  S <- switch(match.arg(sim.mat),
              cts = cts(E = E2, dc = dcCTS),
              srs = srs(E = E2, dc = dcSRS, R = R),
              asrs = asrs(E = E2, dc = dcSRS))
  LCE_cl <- consensus_class(S, k)
  return(LCE_cl)
}
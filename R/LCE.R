#' Linkage Clustering Ensemble
#' 
#' Generate a cluster assignment from a CTS, SRS, or ASRS similarity matrix.
#' 
#' @param E is an array of clustering results. An error is thrown if there are 
#'   missing values. \code{\link{imputeMissing}} can be used beforehand.
#' @param data original data matrix with rows as samples, columns as variables
#' @param k requested number of clusters
#' @param dcCTS decay constant for CTS matrix
#' @param dcSRS decay constant for SRS matrix
#' @param dcASRS decay constant for ASRS matrix
#' @param R number of repetitions for SRS matrix
#' @param sim.mat similarity matrix; choices are "cts", "srs", "asrs".
#' @author Johnson Liu
#' @return a vector containing the cluster assignment from either the CTS, SRS, 
#'   or ASRS similarity matrices
#' @export
#' @examples
#' data(hgsc)
#' dat <- t(hgsc[, -1])[1:100, 1:50]
#' x <- ConClust(dat, nc = 4, reps = 4,
#'               method = c("nmfEucl", "hcAEucl", "hcDianaEucl"), save = FALSE)
#' \dontrun{
#' LCE(E = x, data = dat, k = 4, sim.mat = "asrs")
#' }
#' 
#' x_imputed <- imputeMissing(x, dat, imputeALL = TRUE)$E_imputed2
#' LCE(E = x_imputed, data = dat, k = 4, sim.mat = "cts")
LCE <- function(E, data, k, dcCTS = 0.8, dcSRS = 0.8, dcASRS = 0.8, R = 10,
                sim.mat = c("cts", "srs", "asrs")) {
  assertthat::assert_that(is.array(E),
                          dcCTS >= 0 && dcCTS <= 1,
                          dcSRS >= 0 && dcSRS <= 1,
                          dcASRS >= 0 && dcASRS <= 1)
  # Check that the Cluster matrix is complete otherwise return Error
  if (anyNA(E)) stop("'E' must be complete for LCE algorithm.")
  S <- switch(match.arg(sim.mat,c("cts","asrs","srs")),
              cts = cts(E = E, dc = dcCTS),
              srs = srs(E = E, dc = dcSRS, R = R),
              asrs = asrs(E = E, dc = dcSRS))
  LCE_cl <- consensus_class(S, k)
  return(LCE_cl)
}
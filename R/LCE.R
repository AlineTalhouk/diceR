#' Linkage Clustering Ensemble
#' 
#' Generate CTS, SRS, and ASRS similarity matrices.
#' 
#' @param E is an array of clustering results. Missing entries are imputed using
#'   \code{imputeMissing}
#' @param data original data matrix with rows as samples, columns as variables
#' @param dcCTS decay constant for CTS matrix
#' @param dcSRS decay constant for SRS matrix
#' @param dcASRS decay constant for ASRS matrix
#' @param R number of repetitions for SRS matrix
#' @param sim.mat vector of similarity matrices to compute. Choices are "cts",
#'   "srs", "asrs".
#' @param is.relabelled logical; defaults to \code{TRUE}, but if \code{FALSE}
#'   the data will be relabelled
#' @author Johnson Liu
#' @return a list containing the CTS, SRS, and ASRS matrices (as specified)
#' @export
#' @examples
#' data(hgsc)
#' dat <- t(hgsc[, -1])[1:100, 1:50]
#' x <- ConClust(dat, nc = 4, reps = 4,
#'               method = c("nmfEucl", "hcAEucl", "hcDianaEucl"), save = FALSE)
<<<<<<< HEAD
#' y <- LCE(E = x, data = dat, dcCTS = 0.8, dcSRS = 0.8, dcASRS = 0.8, R = 10,
#' sim.mat = "asrs")
LCE <- function(E, data, dcCTS = 0.8, dcSRS = 0.8, dcASRS = 0.8, R = 10,
                sim.mat = c("cts", "srs", "asrs"), is.relabelled = TRUE) {
  assertthat::assert_that(is.array(E), dcCTS >= 0 && dcCTS <= 1,
                          dcASRS >= 0 && dcASRS <= 1, dcSRS >= 0 && dcSRS <= 1,
                          is_pos_int(R))
=======
#' y <- link_clust(E = x, dcCTS = 0.8, dcSRS = 0.8, dcASRS = 0.8, R = 10)
LCE <- function(E2, k, dcCTS = 0.8, dcSRS = 0.8, dcASRS = 0.8, R = 10,
                sim.mat ="cts") {
>>>>>>> master
  sm.choices <- c("cts", "srs", "asrs")
  if (!all(sim.mat %in% sm.choices)) {
    stop("At least one of 'sim.mat' is not 'cts', 'srs', or 'asrs'.")  
  } else {
    sm <- match.arg(sim.mat, sm.choices, several.ok = TRUE)
  }
<<<<<<< HEAD
  E <- imputeMissing(E, data, imputeALL = TRUE)[, , 1]
  if ("cts" %in% sm)
    CTS <- cts(E = E, dc = dcCTS)
  else
    CTS <- NULL
  if ("srs" %in% sm)
    SRS <- srs(E = E, dc = dcSRS, R = R)
  else
    SRS <- NULL
  if ("asrs" %in% sm)
    ASRS <- asrs(E = E, dc = dcASRS)
  else
    ASRS <- NULL
  out <- list(CTS = CTS, SRS = SRS, ASRS = ASRS)
  return(Filter(Negate(is.null), out))
}
=======
  assertthat::assert_that(is.matrix(E2), dcCTS >= 0 && dcCTS <= 1,
                          dcASRS >= 0 && dcASRS <= 1, dcSRS >= 0 && dcSRS <= 1)
  
    S <- switch (sm,
                 cts = cts(E=E2,dc=dcCTS),
                 asrs = asrs(E=E2,dc=dcSRS),
                 src = srs(E=E2,dc=dcSRS,R=R)
    )
  
  LCE_cl <- consensus_class(S,k)
  return(LCE_cl)
}
>>>>>>> master

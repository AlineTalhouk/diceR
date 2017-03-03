#' Simulate and select null distributions on empirical gene-gene correlations
#'
#' Using a principal component constructed from the sample space, we simulate
#' null distributions with univariate Normal distributions using
#' \code{pcn_simulate}. Then a subset of these distributions is chosen using
#' \code{pcn_select}.
#'
#' @param data data matrix with rows as samples, columns as features
#' @param n.sim The number of simulated datasets to simulate
#'
#' @return \code{pcn_simulate} returns a list of length \code{n.sim}. Each
#'   element is a simulated matrix using this "Principal Component Normal" (pcn)
#'   procedure.
#' @author Derek Chiu
#' @name pcn
#' @export
#' @examples
#' set.seed(9)
#' A <- matrix(rnorm(300), nrow = 20)
#' pc.dat <- pcn_simulate(A, n.sim = 50)
#' cl <- sample(1:4, 20, replace = TRUE)
#' pc.select <- pcn_select(pc.dat, cl, "rep")
pcn_simulate <- function(data, n.sim = 50) {
  pc <- stats::princomp(data)
  Yns <- replicate(n.sim, as.matrix(
    purrr::map_df(pc$sdev, ~ stats::rnorm(pc$n.obs, sd = .x))),
    simplify = FALSE)
  Qns <- purrr::map(Yns, ~ .x %*% t(pc$loadings))
  return(Qns)
}

#' @param data.sim an object from \code{pcn_simulate}
#' @param cl vector of cluster memberships
#' @param type select either the representative dataset ("rep") or a range of
#'   datasets ("range")
#' @param int every \code{int} data sets from median-ranked \code{data.sim} are
#'   taken. Defaults to 5.
#' @return \code{pcn_select} returns a list with elements
#' \item{ranks}{When \code{type = "range"}, ranks of each extracted dataset
#' shown}
#' \item{ind}{index of representative simulation}
#' \item{dat}{simulation data representation of all in pcNormal}
#' @rdname pcn
#' @export
pcn_select <- function(data.sim, cl, type = c("rep", "range"), int = 5) {
  ss <- purrr::map_df(data.sim, sil_widths, cl = cl)
  dists <- apply(ss, 1, function(x)
    stats::dist(rbind(c(stats::median(ss$fN, na.rm = TRUE),
                        stats::median(ss$aP, na.rm = TRUE)),
                      x)))
  type <- match.arg(type)
  switch(type,
         rep = {
           ind <- which.min(dists)
           sim <- data.sim[[ind]]
           return(list(ind = ind, dat = sim))
         },
         range = {
           rks <- seq(1 + int, length(data.sim), int)
           ind <- order(dists)[rks]
           sim <- data.sim[ind]
           return(list(ranks = rks, ind = ind, dat = sim))
         })
}

#' Computes the fraction of negative silhouette widths and average of positive
#' silhouette widths.
#'
#' @param data data matrix
#' @param cl integer vector of cluster memberships
#' @noRd
sil_widths <- function(data, cl) {
  sw <- cluster::silhouette(cl, stats::dist(data))[, "sil_width"]
  return(data.frame(fN = mean(sw < 0), 
                    aP = mean(sw[sw > 0])))
}
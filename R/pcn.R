#' Simulate and select null distributions on empirical gene-gene correlations
#'
#' Using a principal component constructed from the sample space, we simulate
#' null distributions with univariate Normal distributions using `pcn_simulate`.
#' Then a subset of these distributions is chosen using `pcn_select`.
#'
#' @param data data matrix with rows as samples, columns as features
#' @param n.sim The number of simulated datasets to simulate
#'
#' @return `pcn_simulate` returns a list of length `n.sim`. Each element is a
#'   simulated matrix using this "Principal Component Normal" (pcn) procedure.
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
  Yns <- purrr::rerun(n.sim, as.matrix(
    purrr::map_df(pc$sdev, ~ stats::rnorm(pc$n.obs, sd = .x))
  ))
  Qns <- purrr::map(Yns, ~ .x %*% t(pc$loadings))
  return(Qns)
}

#' @param data.sim an object from `pcn_simulate`
#' @param cl vector of cluster memberships
#' @param type select either the representative dataset ("rep") or a range of
#'   datasets ("range")
#' @param int every `int` data sets from median-ranked `data.sim` are taken.
#'   Defaults to 5.
#' @return `pcn_select` returns a list with elements
#' * `ranks`: When `type = "range"`, ranks of each extracted dataset shown
#' * `ind`: index of representative simulation
#' * `dat`: simulation data representation of all in pcNormal
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
           dat <- data.sim[[ind]]
           dplyr::lst(ind, dat)
         },
         range = {
           rks <- seq(1 + int, length(data.sim), int)
           ind <- order(dists)[rks]
           dat <- data.sim[ind]
           dplyr::lst(rks, ind, dat)
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
  data.frame(fN = mean(sw < 0),
             aP = mean(sw[sw > 0]))
}

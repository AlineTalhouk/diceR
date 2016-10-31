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
  Yns <- replicate(n.sim, sapply(pc$sdev, function(x)
    stats::rnorm(pc$n.obs, sd = x)),
    simplify = FALSE)
  Qns <- lapply(Yns, function(x) x %*% t(pc$loadings))
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
  ss <- plyr::ldply(data.sim, sil_widths, cl)
  type <- match.arg(type)
  switch(type,
         rep = {
           ind <- which.min(apply(ss, 1, function(x)
             stats::dist(rbind(c(stats::median(ss$fN), stats::median(ss$aP)),
                               x))))
           sim <- data.sim[[ind]]
           return(list(ind = ind, dat = sim))
         },
         range = {
           rks <- seq(1 + int, length(data.sim), int)
           ord <- order(apply(ss, 1, function(x)
             stats::dist(rbind(c(stats::median(ss$fN), stats::median(ss$aP)),
                               x))))
           ind <- ord[rks]
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
  s <- cluster::silhouette(cl, stats::dist(data))
  fN <- sum(s[, "sil_width"] < 0) / nrow(s)
  aP <- mean(s[, "sil_width"][s[, "sil_width"] > 0])
  return(data.frame(fN, aP))
}

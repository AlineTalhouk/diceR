#' Significant Testing of Clustering Results
#'
#' Uses the SigClust K-Means algorithm to assess significance of clustering
#' results.
#'
#' This function is a wrapper for the original \code{\link[sigclust]{sigclust}},
#' except that an additional parameter \code{k} is allows testing against any
#' number of clusters. In addition, the default type of covariance estimation is
#' also different.
#'
#' @param x data matrix, samples are rows and features are columns
#' @param k cluster size to test against
#' @param nsim number of simulations
#' @param nrep See \code{\link[sigclust]{sigclust}} for details.
#' @param labflag See \code{\link[sigclust]{sigclust}} for details.
#' @param label true class label. See \code{\link[sigclust]{sigclust}} for
#'   details.
#' @param icovest type of covariance matrix estimation
#'
#' @return An object of class \code{sigclust}. See
#'   \code{\link[sigclust]{sigclust}} for details.
#' @author Hanwen Huang: \email{hanwenh@email.unc.edu}; Yufeng Liu:
#'   \email{yfliu@email.unc.edu}; J. S. Marron: \email{marron@email.unc.edu}
#' @references Liu, Yufeng, Hayes, David Neil, Nobel, Andrew and Marron, J. S,
#'   2008, \emph{Statistical Significance of Clustering for High-Dimension,
#'   Low-Sample Size Data}, \emph{Journal of the American Statistical
#'   Association} \bold{103}(483) 1281--1293.
#' @export
#'
#' @examples
#' data(hgsc)
#' dat <- hgsc[1:100, 1:50]
#' nk <- 4
#' cc <- consensus_cluster(dat, nk = nk, reps = 5, algorithms = "pam",
#' progress = FALSE)
#' cl.mat <- consensus_combine(cc, element = "class")
#' lab <- cl.mat$`4`[, 1]
#' set.seed(1)
#' str(sigclust(x = dat, k = nk, nsim = 50, labflag = 1, label = lab))
sigclust <- function(x, k, nsim, nrep = 1, labflag = 0, label = 0,
                     icovest = 2) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  if (n > 1) {
    x <- as.matrix(x)
    if (labflag == 0) {
      xclust <- .cluster(x, k)
      for (i in 1:nrep) {
        clust.temp <- .cluster(x, k)
        if (clust.temp$cindex < xclust$cindex)
          xclust <- clust.temp
        xcindex <- xclust$cindex
      }
    }
    xvareigen <- sigclust::.vareigen(x, n, p, icovest)
    simcindex <- rep(0, nsim)
    for (i in 1:nsim) {
      xsim <- .simnull(xvareigen$vsimeigval, n, p, k)
      simcindex[i] <- xsim$cindex
    }
    if (labflag == 0) {
      index <- simcindex <= xclust$cindex
      mindex <- mean(simcindex)
      sindex <- stats::sd(simcindex)
      pval <- sum(index) / nsim
      pvalnorm <- stats::pnorm(xclust$cindex, mindex, sindex)
    }
    if (labflag == 1) {
      meanpl <- vapply(sort(unique(label)),
                       function(y) colMeans(x[label == y, , drop = FALSE]),
                       FUN.VALUE = double(p))
      txdiffl <- lapply(sort(unique(label)),
                        function(y) t(x[label == y, ]) - meanpl[, y])
      withinsum <- sum(vapply(txdiffl,
                              function(y) sum(y ^ 2),
                              FUN.VALUE = double(1)))
      meanp <- colMeans(x)
      tx <- t(x)
      txdiff <- tx - meanp
      totalsum <- sum(txdiff ^ 2)
      cindexlab <- withinsum / totalsum
      index <- simcindex <= cindexlab
      mindex <- mean(simcindex)
      sindex <- stats::sd(simcindex)
      pval <- sum(index) / nsim
      pvalnorm <- stats::pnorm(cindexlab, mindex, sindex)
      xcindex <- cindexlab
    }
    return(methods::new("sigclust", raw.data = x, veigval = xvareigen$veigval,
                        vsimeigval = xvareigen$vsimeigval,
                        simbackvar = xvareigen$simbackvar,
                        icovest = icovest, nsim = nsim, simcindex = simcindex,
                        pval = pval, pvalnorm = pvalnorm, xcindex = xcindex))
  }
  else {
    print("Only one sample left, no need for clustering!")
    return(0)
  }
}

#' @noRd
.cluster <- function(x, k) {
  clust <- stats::kmeans(x, k)
  withinsum <- sum(clust$withinss)
  meanp <- colMeans(x)
  tx <- t(x)
  txdiff <- tx - meanp
  totalsum <- sum(txdiff ^ 2)
  cindex <- withinsum / totalsum
  list(clust = clust, cindex = cindex)
}

#' @noRd
.simnull <- function(vsimeigval, n, p, k) {
  simnorm <- t(replicate(n, stats::rnorm(p, sd = sqrt(vsimeigval))))
  simclust <- .cluster(simnorm, k)
  list(cindex = simclust$cindex)
}

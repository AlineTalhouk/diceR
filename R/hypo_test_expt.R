# Test hypothesis of unimodal distribution vs. k clusters
# library(diceR)
# library(sigclust)
# data(hgsc)
# dat <- t(hgsc[, -1])
# set.seed(1)
# nk <- 4
# label <- sample(seq_len(nk), nrow(dat), replace = TRUE)
# cc <- consensus_cluster(dat, nk = 4, algorithms = c("hc", "diana", "km", "pam"))
# cl.mat <- consensus_combine(cc, element = "class")
# pvalue <- my_sigclust(dat, k = nk, 1000, labflag = 1, label = cl.mat$`4`[, 4])
# plot(pvalue)

# Modified functions to account for k > 2
my_sigclust <- function(x, k, nsim, nrep = 1, labflag = 0, label = 0, icovest = 1) {
  n <- dim(x)[1]
  p <- dim(x)[2]
  if (n > 1) {
    x <- as.matrix(x)
    if (labflag == 0) {
      xclust <- my_.cluster(x, n, p, k)
      for (i in 1:nrep) {
        clust.temp <- my_.cluster(x, n, p, k)
        if (clust.temp$cindex < xclust$cindex) 
          xclust <- clust.temp
        xcindex <- xclust$cindex
      }
    }
    xvareigen <- sigclust::.vareigen(x, n, p, icovest)
    simcindex <- rep(0, nsim)
    for (i in 1:nsim) {
      xsim <- my_.simnull(xvareigen$vsimeigval, n, p, k)
      simcindex[i] <- xsim$cindex
    }
    if (labflag == 0) {
      index <- (simcindex <= xclust$cindex)
      mindex <- mean(simcindex)
      sindex <- sd(simcindex)
      pval <- sum(index)/nsim
      pvalnorm <- pnorm(xclust$cindex, mindex, sindex)
    }
    if (labflag == 1) {
      meanpl <- sapply(sort(unique(label)), function(y) colMeans(x[label == y, ]))
      txdiffl <- sapply(sort(unique(label)), function(y) t(dat[label == y, ]) - meanpl[, y])
      withinsum <- sum(sapply(txdiffl, function(y) sum(y^2)))
      # meanp1 <- colMeans(x[label == 1, ])
      # txdiff1 <- t(x[label == 1, ]) - meanp1
      # meanp2 <- colMeans(x[label == 2, ])
      # txdiff2 <- t(x[label == 2, ]) - meanp2
      # withinsum <- sum(txdiff1^2) + sum(txdiff2^2)
      meanp <- colMeans(x)
      tx <- t(x)
      txdiff <- tx - meanp
      totalsum <- sum(txdiff^2)
      cindexlab <- withinsum/totalsum
      index <- (simcindex <= cindexlab)
      mindex <- mean(simcindex)
      sindex <- sd(simcindex)
      pval <- sum(index)/nsim
      pvalnorm <- pnorm(cindexlab, mindex, sindex)
      xcindex <- cindexlab
    }
    return(new("sigclust", raw.data = x, veigval = xvareigen$veigval, 
               vsimeigval = xvareigen$vsimeigval, simbackvar = xvareigen$simbackvar, 
               icovest = icovest, nsim = nsim, simcindex = simcindex, 
               pval = pval, pvalnorm = pvalnorm, xcindex = xcindex))
  }
  else {
    print("Only one sample left, no need for clustering!")
    return(0)
  }
}

my_.cluster <- function(x, n, p, k) {
  if (n > 1) {
    x <- as.matrix(x)
    if (dim(x)[1] == n & dim(x)[2] == p) {
      clust <- kmeans(x, k)
      withinsum <- sum(clust$withinss)
      meanp <- colMeans(x)
      tx <- t(x)
      txdiff <- tx - meanp
      totalsum <- sum(txdiff^2)
      cindex <- withinsum / totalsum
      list(clust = clust, cindex = cindex)
    }
    else {
      print("Wrong size of matrix x!")
      return(0)
    }
  }
  else {
    print("Only one sample left, no need for clustering!")
    return(0)
  }
}

my_.simnull <- function(vsimeigval, n, p, k) {
  simnorm <- matrix(0, n, p)
  for (i in 1:n) {
    simnorm[i, ] <- rnorm(p, sd = sqrt(vsimeigval))
  }
  simclust <- my_.cluster(simnorm, n, p, k)
  list(cindex = simclust$cindex)
}
library(diceR)
data(hgsc)
dat <- t(hgsc[,-1])
s <- 1
algs <- "nmfbrunet"
pr <- "cs"
k = 4
r = 1
print(sessionInfo())  # added this line
ssclust <- switch(
  algs,
  nmfbrunet = consensus_cluster(
    data = dat,
    nk = k,
    reps = r,
    algorithms = "nmf",
    nmf.method = "brunet",
    prep.data = "none",
    seed.data = s,
    save = FALSE
  )
)

context("sigclust")

data(hgsc)
dat <- t(hgsc[, -1])[1:100, 1:50]
nk <- 4
cc <- consensus_cluster(dat, nk = 4, reps = 2, algorithms = "hc",
                        progress = FALSE)
cl.mat <- consensus_combine(cc, element = "class")

test_that("sigclust works regardless of labflag param", {
  set.seed(1)
  expect_error(sigclust(x = dat, k = nk, nsim = 1000, labflag = 0,
                        label = cl.mat$`4`[, 1], icovest = 2), NA)
  expect_error(sigclust(x = dat, k = nk, nsim = 1000, labflag = 1,
                        label = cl.mat$`4`[, 1], icovest = 2), NA)
})

test_that("cindex loop and condition work", {
  set.seed(4)
  expect_error(sigclust(x = dat, k = nk, nsim = 1000, nrep = 100, labflag = 0,
                        label = cl.mat$`4`[, 1], icovest = 2), NA)
})

test_that("one sample matrix isn't validated", {
  expect_output(sigclust(x = dat[1, , drop = FALSE], k = nk, nsim = 1000))
})
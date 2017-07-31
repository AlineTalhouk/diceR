context("Signifance Testing")

data(hgsc)
dat <- hgsc[1:100, 1:50]
nk <- 4
cc <- consensus_cluster(dat, nk = 4, reps = 2, algorithms = "hc",
                        progress = FALSE)
cl.mat <- consensus_combine(cc, element = "class")
ref_label <- cl.mat$`4`[, 1]

test_that("sigclust works regardless of labflag param", {
  set.seed(1)
  expect_error(sigclust(x = dat, k = nk, nsim = 10, labflag = 0,
                        label = ref_label), NA)
  expect_error(sigclust(x = dat, k = nk, nsim = 10, labflag = 1,
                        label = ref_label), NA)
})

test_that("cindex loop and condition work", {
  set.seed(4)
  expect_error(sigclust(x = dat, k = nk, nsim = 10, nrep = 100, labflag = 0,
                        label = ref_label), NA)
})

test_that("one sample matrix isn't validated", {
  expect_output(sigclust(x = dat[1, , drop = FALSE], k = nk, nsim = 10))
})

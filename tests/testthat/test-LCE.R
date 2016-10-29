
context("Linkage Cluster Ensemble")

test_that("Check LCE with hgsc data with 3 ConClust algorithms", {
  data(hgsc)
  dat <- t(hgsc[, -1])[1:200, 1:100]
  x <- ConClust(dat, nc = 4, reps = 4, progress = FALSE,
                method = c("nmfEucl", "hcAEucl", "hcDianaEucl"))
  y_cts <- LCE(E = x, data = dat, k = 4, dcCTS = 0.8, dcSRS = 0.8,
           dcASRS = 0.8, R = 5, sim.mat = "cts")
  y_srs <- LCE(E = x, data = dat, k = 4, dcCTS = 0.8, dcSRS = 0.8,
               dcASRS = 0.8, R = 5, sim.mat = "srs")
  y_asrs <- LCE(E = x, data = dat, k = 4, dcCTS = 0.8, dcSRS = 0.8,
               dcASRS = 0.8, R = 5, sim.mat = "asrs")
  expect_length(y_cts, 200)
  expect_length(y_srs, 200)
  expect_length(y_asrs, 200)
  expect_is(y_cts, "factor")
  expect_is(y_srs, "factor")
  expect_is(y_asrs, "factor")
})
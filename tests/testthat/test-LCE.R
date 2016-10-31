
context("Linkage Cluster Ensemble")

test_that("Check LCE with hgsc data with 3 ConClust algorithms", {
  data(hgsc)
  dat <- t(hgsc[, -1])[1:200, 1:100]
  k <- 4
  x <- ConClust(dat, nc = k, reps = 4, progress = FALSE,
                method = c("nmfEucl", "hcAEucl", "hcDianaEucl"))
  x_imputed <- imputeMissing(x, dat)$E_imputed2
  y_cts <- LCE(E = x_imputed, data = dat, k = k, R = 5, sim.mat = "cts")
  y_srs <- LCE(E = x_imputed, data = dat, k = k, R = 5, sim.mat = "srs")
  y_asrs <- LCE(E = x_imputed, data = dat, k = k, R = 5, sim.mat = "asrs")
  expect_length(y_cts, 200)
  expect_length(y_srs, 200)
  expect_length(y_asrs, 200)
  expect_is(y_cts, "factor")
  expect_is(y_srs, "factor")
  expect_is(y_asrs, "factor")
})

context("Consensus functions")

test_that("majority voting works", {
  dt <- array(c(2, 3, 2, 2, 2, 3, 1, 1, 2, 3, 2, 2, 2, 3, 1, 1), c(2, 4, 2))
  expect_equal(majority_voting(dt, is.relabelled = TRUE), c(2, 3))
})

test_that("k modes works with or without missing", {
  x <- array(rep(c(rep(1, 10), rep(2, 10), rep(3, 10)), times = 5), c(30, 6, 5))
  xf <- x
  dim(xf) <- c(30, 30)
  set.seed(1)
  kmo.old <- klaR::kmodes(xf, 3)$cluster
  kmo.new <- k_modes(x)
  expect_equal(kmo.old, kmo.new)
  
  x[3, , 1] <- NA
  kmo.missing <- k_modes(x)
  expect_false(anyNA(kmo.missing))
})

set.seed(1)
x <- replicate(100, rbinom(100, 4, 0.2))
y <- consensus_matrix(x)

test_that("CSPA works", {
  expect_length(CSPA(y, k = 2), nrow(y))
  expect_equal(dplyr::n_distinct(CSPA(y, k = 2)), 2)
  expect_equal(dplyr::n_distinct(CSPA(y, k = 3)), 3)
})

test_that("Check LCE with hgsc data with 3 ConClust algorithms", {
  data(hgsc)
  dat <- t(hgsc[, -1])[1:200, 1:100]
  k <- 4
  x <- ConClust(dat, nk = k, reps = 4, progress = FALSE,
                method = c("nmfEucl", "hcAEucl", "hcDianaEucl"))
  x_imputed <- impute_missing(x, dat)$complete
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
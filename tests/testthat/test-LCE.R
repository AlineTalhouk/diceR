
context("Linkage Cluster Ensemble")

test_that("Check LCE with hgsc data with 3 ConClust algorithms", {
  data(hgsc)
  dat <- t(hgsc[, -1])[1:200, 1:100]
  x <- ConClust(dat, nc = 4, reps = 4, progress = FALSE,
                method = c("nmfEucl", "hcAEucl", "hcDianaEucl"))
  y <- LCE(E = x, data = dat, dcCTS = 0.8, dcSRS = 0.8,
           dcASRS = 0.8, R = 5,
           is.relabelled = FALSE)
  y_cts <- y$CTS
  y_srs <- y$SRS
  y_asrs <- y$ASRS
  expect_equal(dim(y_cts), c(200, 200))
  expect_equal(dim(y_srs), c(200, 200))
  expect_equal(dim(y_asrs), c(200, 200))
  expect_equal(sum(diag(y_cts) != 1), 0)
  expect_equal(sum(diag(y_srs) != 1), 0)
  expect_equal(sum(diag(y_asrs) != 1), 0)
  expect_true(is.numeric(y_cts))
  expect_true(is.numeric(y_srs))
  expect_true(is.numeric(y_asrs))
  expect_equal(sum(!complete.cases(y_cts)), 0)
  expect_equal(sum(!complete.cases(y_srs)), 0)
  expect_equal(sum(!complete.cases(y_asrs)), 0)
})
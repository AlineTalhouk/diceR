
context("Link clustering")

set.seed(1)
E1 <- matrix(rep(sample(1:4, 1000, replace = TRUE)), nrow = 100)
E2 <- matrix(rep(sample(c(1:4, NA), 1000, replace = TRUE)), nrow = 100)

test_that("Check link_clust with 2D complete case matrix", {
  test_lc_2d_complete <- link_clust(E = E1, dcCTS = 0.8, dcSRS = 0.8,
                                    dcASRS = 0.8, R = 10, is.relabelled = FALSE)
  expect_equal(sum(test_lc_2d_complete$cts != cts(E1, 0.8)), 0)
  expect_equal(sum(test_lc_2d_complete$srs != srs(E1, 0.8, 10)), 0)
  expect_equal(sum(test_lc_2d_complete$asrs != asrs(E1, 0.8)), 0)
})

test_that("Check link_clust with 2D incomplete case matrix", {
  test_lc_2d_incomplete <- link_clust(E = E2, dcCTS = 0.8, dcSRS = 0.8,
                                      dcASRS = 0.8, R = 10,
                                      is.relabelled = FALSE)
  expect_true(is.numeric(test_lc_2d_incomplete$cts))
  expect_true(is.numeric(test_lc_2d_incomplete$srs))
  expect_true(is.numeric(test_lc_2d_incomplete$asrs))
  expect_equal(sum(diag(test_lc_2d_incomplete$cts) != 1), 0)
  expect_equal(sum(diag(test_lc_2d_incomplete$srs) != 1), 0)
  expect_equal(sum(diag(test_lc_2d_incomplete$asrs) != 1), 0)
  expect_equal(dim(test_lc_2d_incomplete$cts), c(100, 100))
  expect_equal(dim(test_lc_2d_incomplete$srs), c(100, 100))
  expect_equal(dim(test_lc_2d_incomplete$asrs), c(100, 100))
  expect_equal(sum(!complete.cases(test_lc_2d_incomplete$cts)), 0)
  expect_equal(sum(!complete.cases(test_lc_2d_incomplete$srs)), 0)
  expect_equal(sum(!complete.cases(test_lc_2d_incomplete$asrs)), 0)
})

test_that("Check link_clust with hgsc data with 3 ConClust algorithms", {
  data(hgsc)
  dat <- t(hgsc[, -1])
  x <- ConClust(dat[1:200, 1:100], k = 4, reps = 4, progress = FALSE,
                method = c("nmfEucl", "hcAEucl", "hcDianaEucl"))
  test_lc_3d_incomplete <- link_clust(E = x, dcCTS = 0.8, dcSRS = 0.8,
                                      dcASRS = 0.8, R = 5,
                                      is.relabelled = FALSE)
  expect_equal(dim(test_lc_3d_incomplete$cts), c(200, 200))
  expect_equal(dim(test_lc_3d_incomplete$srs), c(200, 200))
  expect_equal(dim(test_lc_3d_incomplete$asrs), c(200, 200))
  expect_equal(sum(diag(test_lc_3d_incomplete$cts) != 1), 0)
  expect_equal(sum(diag(test_lc_3d_incomplete$asrs) != 1), 0)
  expect_equal(sum(diag(test_lc_3d_incomplete$srs) != 1), 0)
  expect_true(is.numeric(test_lc_3d_incomplete$cts))
  expect_true(is.numeric(test_lc_3d_incomplete$srs))
  expect_true(is.numeric(test_lc_3d_incomplete$asrs))
  expect_equal(sum(!complete.cases(test_lc_3d_incomplete$cts)), 0)
  expect_equal(sum(!complete.cases(test_lc_3d_incomplete$srs)), 0)
  expect_equal(sum(!complete.cases(test_lc_3d_incomplete$asrs)), 0)
})
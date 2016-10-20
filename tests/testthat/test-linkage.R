
context("Linkage")

set.seed(1)
E_LCE <- matrix(rep(sample(1:4, 800, replace = TRUE)), nrow = 100)

test_that("Check linkage with complete linkage", {
  Y <- linkage(stod(asrs(E_LCE, 0.8)), "complete")
  expect_equal(nrow(Y), nrow(E_LCE) - 1)
  expect_equal(abs(sum(Y[, 3])), 37.65919, tolerance = 0.0001)
})

test_that("Check linkage with average linkage", {
  Y <- linkage(stod(asrs(E_LCE, 0.8)), "average")
  expect_equal(nrow(Y), nrow(E_LCE) - 1)
  expect_equal(abs(sum(Y[, 3])), 36.55684, tolerance = 0.0001)
})

test_that("Check linkage with single linkage", {
  Y <- linkage(stod(asrs(E_LCE, 0.8)), "single")
  expect_equal(nrow(Y), nrow(E_LCE) - 1)
  expect_equal(abs(sum(Y[, 3])), 33.87901, tolerance = 0.0001)
})

test_that("Check linkage with wrong inputs", {
  expect_error(linkage(c("a", "b", "c"), "complete"))
  expect_error(linkage(stod(asrs(E_LCE, 0.8)), "bccrc"))
  expect_error(linkage(asrs(E_LCE, 0.8), "average"))
})
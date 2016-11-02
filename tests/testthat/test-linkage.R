
context("Linkage")

set.seed(1)
E_LCE <- matrix(rep(sample(1:4, 800, replace = TRUE)), nrow = 100)

test_that("Check linkage with complete linkage", {
  Y <- linkage(asrs(E_LCE, 0.8), "complete")
  expect_equal(nrow(Y), nrow(E_LCE) - 1)
  expect_equal(abs(sum(Y[, 3])), 37.65919, tolerance = 0.0001)
})

test_that("Check linkage with average linkage", {
  Y <- linkage(asrs(E_LCE, 0.8), "average")
  expect_equal(nrow(Y), nrow(E_LCE) - 1)
  expect_equal(abs(sum(Y[, 3])), 36.55684, tolerance = 0.0001)
})

test_that("Check linkage with single linkage", {
  Y <- linkage(asrs(E_LCE, 0.8), "single")
  expect_equal(nrow(Y), nrow(E_LCE) - 1)
  expect_equal(abs(sum(Y[, 3])), 33.87901, tolerance = 0.0001)
})

test_that("Check linkage with wrong inputs", {
  expect_error(linkage(c("a", "b", "c"), "complete"))
  expect_error(linkage(asrs(E_LCE, 0.8), "bccrc"))
  expect_error(linkage(stod(asrs(E_LCE, 0.8)), "average"))
})


context("stod, distance vector")

test_that("Check stod output", {
  S <- cts(E = E_LCE, dc = 0.8)
  s <- stod(S)
  expect_equal(abs(sum(s)), 1097, tolerance = 0.01)
})

test_that("Check stod with wrong inputs",{
  expect_error(stod(1))
  expect_error(stod(1:2))
  expect_error(stod("vancouver"))
  expect_error(stod(matrix(letters[1:6], nrow = 3)))
})
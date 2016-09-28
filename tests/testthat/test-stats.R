
context("Cluster analysis statistics")

test_that("accuracy works", {
  set.seed(1)
  x <- matrix(rbinom(16, 20, 0.4), nrow = 4)
  expect_error(accuracy(x), NA)
})

test_that("normalized mutual information works", {
  set.seed(1)
  X <- sample(1:4, 100, replace = TRUE)
  Y <- sample(1:4, 100, replace = TRUE)
  expect_error(NMI(X, Y), NA)
})

test_that("PAC can have different bounds", {
  set.seed(1)
  x <- replicate(100, rbinom(100, 4, 0.2))
  y <- consensus_matrix(x)
  expect_error(PAC(y, lower = 0.3, upper = 0.7), NA)
})

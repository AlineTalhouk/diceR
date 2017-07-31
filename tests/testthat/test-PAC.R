context("Proportion of Ambiguous Clusterings")

test_that("PAC can have different bounds", {
  set.seed(1)
  x <- replicate(100, rbinom(100, 4, 0.2))
  y <- consensus_matrix(x)
  expect_error(PAC(y, lower = 0.3, upper = 0.7), NA)
})

test_that("A zero matrix returns NA", {
  zero_mat <- consensus_matrix(matrix(NA, ncol = 3, nrow = 3))
  expect_equal(PAC(zero_mat), NA)
})

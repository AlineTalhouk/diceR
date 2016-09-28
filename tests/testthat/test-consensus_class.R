
context("Consensus class")

set.seed(1)
x <- replicate(100, rbinom(100, 4, 0.2))
y <- consensus_matrix(x)

test_that("consensus class works", {
  expect_length(consensus_class(y, k = 2), nrow(y))
  expect_equal(dplyr::n_distinct(consensus_class(y, k = 2)), 2)
  expect_equal(dplyr::n_distinct(consensus_class(y, k = 3)), 3)
})

test_that("names can be overridden", {
  cc1 <- consensus_class(y, k = 2)
  cc2 <- consensus_class(y, k = 2, names = paste0("S", seq_len(100)))
  expect_true(all(cc1 == cc2))
  expect_null(names(cc1))
  expect_identical(names(cc2),  paste0("S", seq_len(100)))
})

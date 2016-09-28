
context("Hook functions")

library(cluster)
data(votes.repub)

test_that("algorithm hooks work", {
  d <- dist(votes.repub)
  expect_length(agnes_hook(d, k = 2), nrow(votes.repub))
  expect_length(diana_hook(d, k = 2), nrow(votes.repub))
})

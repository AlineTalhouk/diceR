
context("Utility functions")

test_that("data preparation removes rows", {
  set.seed(2)
  x <- replicate(10, rnorm(100, sd = 1))
  expect_lt(nrow(prepare_data(x)), nrow(x))
})

test_that("relabelling outputs a factor", {
  set.seed(2)
  pred <- sample(1:4, 100, replace = TRUE)
  true <- sample(1:4, 100, replace = TRUE)
  expect_is(relabel_class(pred, true), "factor")
})

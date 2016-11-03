
context("External validity indices")

set.seed(1)
E <- matrix(rep(sample(1:4, 1000, replace = TRUE)), nrow = 100)

test_that("normalized mutual information works", {
  set.seed(1)
  x <- sample(1:4, 100, replace = TRUE)
  y <- sample(1:4, 100, replace = TRUE)
  expect_error(ev_nmi(x, y), NA)
})

context("Majority voting and K modes")

test_that("majority voting works", {
  dt <- array(c(2, 3, 2, 2, 2, 3, 1, 1, 2, 3, 2, 2, 2, 3, 1, 1), c(2, 4, 2))
  expect_equal(majority_voting(dt, is.relabelled = TRUE), c(2, 3))
})

test_that("k modes works", {
  x <- array(rep(c(rep(1, 10), rep(2, 10), rep(3, 10)), times = 5), c(30, 6, 5))
  xf <- x
  dim(xf) <- c(30, 30)
  set.seed(1)
  kmo.old <- klaR::kmodes(xf, 3)$cluster
  kmo.new <- k_modes(x, is.relabelled = TRUE, seed = 1)
  expect_equal(kmo.old, kmo.new)
})
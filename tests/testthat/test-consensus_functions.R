
context("Majority voting and K modes")

test_that("majority voting works", {
  dt <- array(c(2, 3, 2, 2, 2, 3, 1, 1, 2, 3, 2, 2, 2, 3, 1, 1), c(2, 4, 2))
  expect_equal(majority_voting(dt, is.relabelled = TRUE), c("2", "3"))
})

test_that("k modes works", {
  set.seed(1)
  x <- rbind(matrix(rbinom(250, 2, 0.25), ncol = 5),
             matrix(rbinom(250, 2, 0.75), ncol = 5))
  colnames(x) <- c("a", "b", "c", "d", "e")
  xp <- x
  dim(xp) <- c(10, 5, 10)
  set.seed(1)
  kmo.old <- klaR::kmodes(xp, 2)$cluster
  set.seed(1)
  kmo.new <- k_modes(xp, is.relabelled = TRUE)
  expect_equal(kmo.old, kmo.new)
})

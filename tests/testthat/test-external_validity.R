
context("External validity indices")

set.seed(1)
E <- matrix(rep(sample(1:4, 1000, replace = TRUE)), nrow = 100)

test_that("Check ev_accuracy works with basic inputs", {
  expect_equal(ev_accuracy(c(1, 2, 1, 1), c(1, 1, 1, 1)), 1)
  expect_equal(ev_accuracy(c(4, 2, 4, 4), c(4, 4, 4, 4)), 1)
  expect_equal(ev_accuracy(c(4, 4, 4, 4), c(4, 2, 4, 4)), 0.75)
  expect_equal(ev_accuracy(c(4, 4, 4, 4, 8), c(4, 2, 4, 4, 6)), 0.8)
})

test_that("Check ev_accuracy throws error with wrong inputs", {
  expect_error(ev_accuracy(c(1, 2, 3, 4)))
  expect_error(ev_accuracy("a"))
  expect_error(ev_accuracy(c("a", "b", "c", "d")))
  expect_equal(ev_accuracy(c("a", "a", "b"), c("a", "b", "c")), 2 / 3)
})

test_that("Check ev_rand with E and in case where AR is non-zero", {
  expect_true(abs(ev_rand(E[, 1], E[, 2])$AR + 0.0017) <= 0.0001)
  expect_true(abs(ev_rand(E[, 1], E[, 2])$HI - 0.2469) <= 0.0001)
  expect_true(abs(ev_rand(E[, 1], E[, 2])$MI - 0.3766) <= 0.0001)
  expect_true(abs(ev_rand(E[, 1], E[, 2])$RI - 0.6234) <= 0.0001)
})

test_that("Check ev_rand with wrong inputs", {
  expect_error(ev_rand(1, 2))
  expect_error(ev_rand(
    data.frame(name = letters[1:4], vals = 1:4),
    data.frame(name = letters[1:4], vals = 1:4)
  ))
  expect_error(ev_rand(c(1, 2, 3, 4),
                       data.frame(name = letters[1:4], vals = 1:4)))
  expect_error(ev_rand(data.frame(name = letters[1:4], vals = 1:4),
                       c(6, 8, 3, 8)))
  expect_error(ev_rand(E[1:10, 1], E[, 2]))
  expect_error(ev_rand(matrix(c(1, 2, 3), ncol = 3),
                       matrix(c(1, 2, 3), ncol = 3)))
})


test_that("normalized mutual information works", {
  set.seed(1)
  x <- sample(1:4, 100, replace = TRUE)
  y <- sample(1:4, 100, replace = TRUE)
  expect_error(ev_nmi(x, y), NA)
})

test_that("accuracy works", {
  set.seed(1)
  x <- matrix(rbinom(16, 20, 0.4), nrow = 4)
  expect_error(accuracy(x), NA)
})
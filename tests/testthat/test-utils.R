context("Utility functions")

cl <- seq_len(4)

test_that("relabelling outputs a integer", {
  set.seed(2)
  pred <- sample(cl, 100, replace = TRUE)
  true <- sample(cl, 100, replace = TRUE)
  expect_is(relabel_class(pred, true), "integer")
})

test_that("flatten uses first clustering as reference if not relabelled", {
  E <- matrix(sample(cl, 500, replace = TRUE), ncol = 5)
  E4 <- array(sample(cl, 2000, replace = TRUE), dim = c(100, 5, 2, 2))
  expect_error(flatten_E(E, is.relabelled = FALSE), NA)
  expect_error(flatten_E(E4, is.relabelled = FALSE), NA)
})

test_that("Check is_pos_int works", {
  expect_true(is_pos_int(3))
  expect_false(is_pos_int(-3))
  expect_true(is_pos_int(1e6))
  expect_false(is_pos_int(3.6))
  expect_false(is_pos_int(3.21e1))
  expect_false(is_pos_int(-3.7))
  expect_false(is_pos_int(-3.21e1))
  expect_false(is_pos_int(0))
})

test_that("Error if input is not a number or more than one element", {
  expect_error(is_pos_int(c(1, 2, 3)))
  expect_error(is_pos_int("a"))
  expect_error(is_pos_int(matrix(c(1, 2, 3), nrow = 1)))
})

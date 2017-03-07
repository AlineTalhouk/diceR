
context("Utility functions")

test_that("relabelling outputs a integer", {
  set.seed(2)
  pred <- sample(1:4, 100, replace = TRUE)
  true <- sample(1:4, 100, replace = TRUE)
  expect_is(relabel_class(pred, true), "integer")
})

test_that("flatten uses first clustering as reference if not relabelled", {
  E <- matrix(sample(1:4, 500, replace = TRUE), ncol = 5)
  E4 <- array(sample(1:4, 2000, replace = TRUE), dim = c(100, 5, 2, 2))
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

test_that("Check sortMatrixRowWise", {
  mat <- matrix(c(8, 63, 93, 7, 30, 16, 5, 22, 48, 32, 89, 56, 1, 83, 76, 9),
                ncol = 4)
  expect_equal(sum(!sortMatrixRowWise(mat, "ascending") == 
                     matrix(c(1, 16, 5, 7, 8, 32, 76, 9, 30, 63,
                              89, 22, 48, 83, 93, 56), ncol = 4)), 0)
  expect_equal(sum(!sortMatrixRowWise(mat, "descending") ==
                     matrix(c(48, 83, 93, 56, 30, 63, 89, 22, 8,
                              32, 76, 9, 1, 16, 5, 7), ncol = 4)), 0)
  mat_reps <- matrix(c(8, 63, 93, 7, 30, 32, 5, 22, 48,
                       32, 89, 56, 1, 32, 76, 9), ncol = 4)
  expect_equal(sum(
    !sortMatrixRowWise(mat_reps, "ascending") ==
      matrix(c(1, 32, 5, 7, 8, 32, 76, 9, 30, 32, 89, 22, 48, 63, 93, 56),
             ncol = 4)), 0)
  expect_equal(sum(
    !sortMatrixRowWise(mat_reps, "descending") ==
      matrix(c(1, 32, 5, 7, 8, 32, 76, 9, 30, 32, 89, 22, 48, 63, 93, 56),
             ncol = 4)), 14)
})

context("Utility functions")

test_that("relabelling outputs a factor", {
  set.seed(2)
  pred <- sample(1:4, 100, replace = TRUE)
  true <- sample(1:4, 100, replace = TRUE)
  expect_is(relabel_class(pred, true), "factor")
})

test_that("flatten uses first clustering as reference if not relabelled", {
  E <- matrix(sample(1:4, 500, replace = TRUE), ncol = 5)
  E4 <- array(sample(1:4, 2000, replace = TRUE), dim = c(100, 5, 2, 2))
  expect_error(flatten_E(E, is.relabelled = FALSE), NA)
  expect_error(flatten_E(E4, is.relabelled = FALSE), NA)
})

a <- matrix(c(60, 17, 58, 62, 81, 11, 32, 7, 28, 85, 80, 15, 19, 50, 45, 40,
              88, 31, 84, 30, 99, 94, 61, 55, 27), ncol = 5)

test_that("Check coord works with a matrix without duplicated entry", {
  expect_true(coord(a, 60)$rows == 1)
  expect_true(coord(a, 60)$cols == 1)
  expect_true(coord(a, 85)$rows == 5)
  expect_true(coord(a, 85)$cols == 2)
  expect_true(coord(a, 84)$rows == 4)
  expect_true(coord(a, 84)$cols == 4)
  expect_true(coord(a, 94)$rows == 2)
  expect_true(coord(a, 94)$cols == 5)
})

a_wd <- matrix(c(60, 17, 58, 62, 81, 11, 32, 7, 28, 85, 80, 15, 19, 50, 45, 40,
                 88, 32, 84, 30, 32, 94, 61, 85, 27), ncol = 5)

test_that("Check getRowColNumbers works with matrix with duplicated values", {
  expect_equal(coord(a_wd, 32)$rows, c(2, 3, 1))
  expect_equal(coord(a_wd, 32)$cols, c(2, 4, 5))
  expect_equal(coord(a_wd, 85)$rows, c(5, 4))
  expect_equal(coord(a_wd, 85)$cols, c(2, 5))
})

test_that("Check getRowColNumbers throw error with wrong type of input", {
  expect_error(coord(a, 101))
  expect_error(coord(a, c(32, 80)))
  expect_error(coord(a, "101"))
  expect_error(coord(data.frame(val1 = c(1, 2, 3), val2 = c(4, 5, 6)), 1))
  expect_error(coord(matrix(letters[1:25], ncol = 5), "a"))
  expect_error(coord(matrix(letters[1:25], ncol = 5), 1))
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
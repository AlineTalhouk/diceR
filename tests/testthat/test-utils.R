
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

a <- matrix(c(60, 17, 58, 62, 81, 11, 32, 7, 28, 85, 80, 15, 19, 50, 45, 40,
              88, 31, 84, 30, 99, 94, 61, 55, 27), ncol = 5)

test_that("Check colMax works with a matrix", {
  expect_equal(colMax(matrix(1:16, ncol = 4)), c(4, 8, 12, 16))
  expect_equal(colMax(a), c(81, 85, 80, 88, 99))
})

test_that("Check colMax throws error with wrong type of input", {
  expect_error(colMax(data.frame(
    names = c("a", "b", "c"), vals = c(17, 18, 19)
  )))
  expect_error(colMax(data.frame(
    states = c("CALI", "TEXS", "NEWY"),
    cities = c("san fran", "houston", "new york")
  )))
  expect_error(colMax(data.frame(
    val1 = c(1, 2, 3), val2 = c(3, 6, 9)
  )))
})

test_that("Check colMin works with a matrix", {
  expect_equal(colMin(matrix(1:16, ncol = 4)), c(1, 5, 9, 13))
  expect_equal(colMin(a), c(17, 7, 15, 30, 27))
})

test_that("Check colMin throws error with wrong type of input",{
  expect_error(colMin(data.frame(
    names = c("a", "b", "c"), vals = c(17, 18, 19)
  )))
  expect_error(colMin(data.frame(
    states = c("CALI", "TEXS", "NEWY"),
    cities = c("san fran", "houston", "new york")
  )))
  expect_error(colMin(data.frame(
    val1 = c(1, 2, 3), val2 = c(4, 5, 6)
  )))
})

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

a_wd <- matrix(c(60, 17, 58, 62, 81, 11, 32, 7, 28, 85, 80, 15, 19, 50, 45,
                 40, 88, 32, 84, 30, 32, 94, 61, 85, 27), ncol = 5)

test_that("Check getRowColNumbers works with matrix with duplicated values",{
  expect_equal(coord(a_wd, 32)$rows, c(2, 3, 1))
  expect_equal(coord(a_wd, 32)$cols, c(2, 4, 5))
  expect_equal(coord(a_wd, 85)$rows, c(5, 4))
  expect_equal(coord(a_wd, 85)$cols, c(2, 5))
})

test_that("Check getRowColNumbers throw error with wrong type of input",{
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

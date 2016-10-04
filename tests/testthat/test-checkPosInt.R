
context("Check positive integer")

test_that("Check checkPosInt works", {
  expect_true(checkPosInt(3))
  expect_false(checkPosInt(-3))
  expect_true(checkPosInt(1e6))
  expect_false(checkPosInt(3.6))
  expect_false(checkPosInt(3.21e1))
  expect_false(checkPosInt(-3.7))
  expect_false(checkPosInt(-3.21e1))
})

test_that("Check if input is not a number gives an error", {
  expect_error(checkPosInt(c(1, 2, 3)))
  expect_error(checkPosInt("a"))
  expect_error(checkPosInt(matrix(c(1, 2, 3), nrow = 1)))
})

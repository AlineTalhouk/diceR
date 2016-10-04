context("Column maximums")

test_that("Check colMax works with a matrix", {

  expect_equal(colMax(matrix(1:16, ncol = 4)), c(4, 8, 12, 16))
  expect_equal(colMax(matrix(
    c(
      60,
      17,
      58,
      62,
      81,
      11,
      32,
      7,
      28,
      85,
      80,
      15,
      19,
      50,
      45,
      40,
      88,
      31,
      84,
      30,
      99,
      94,
      61,
      55,
      27
    ),
    ncol = 5
  )),
  c(81, 85, 80, 88, 99))
})

test_that("Check colMax throws error with wrong type of input",{
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

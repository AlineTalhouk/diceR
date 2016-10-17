
context("stod, distance vector")

set.seed(1)
E <-
  matrix(rep(sample(1:4, 1000, replace = TRUE)), nrow = 100, byrow = FALSE)
dc <- 0.8

test_that("Check stod", {
  S <- cts(E=E, dc = dc)
  s <- stod(S)
  expect_true(abs(sum(s) - 1115) < 1)

})

test_that("Check stod with wrong inpts",{
  expect_error(stod(1))
  expect_error(stod(c(1, 2)))
  expect_error(stod("vancouver"))
  expect_error(stod(matrix(c(
    "a", "b", "c", "d", "e", "f"
  ), nrow = 3)))
})

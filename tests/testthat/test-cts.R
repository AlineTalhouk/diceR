
context("Connected triple based similarity matrix")

data("E_LCE")

test_that("Check cts works", {

  CTS <- cts(E_LCE, 0.8)
  expect_equal(sum(!diag(CTS) == 1), 0)
  expect_true(abs(sum(CTS) - 4076.2) < 0.1)

})

test_that("Check cts throws error with wrong type of input",{
  expect_error(cts(E_LCE))
  expect_error(cts(E_LCE, 1.2))
  expect_error(cts(E_LCE, -9))
  expect_error(cts(c(1, 1, 2, 3), 0.8))
})

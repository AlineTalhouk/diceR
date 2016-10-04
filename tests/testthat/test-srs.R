
context("SRS, simrank based similarity matrix")

data("E_LCE")

test_that("Check srs works", {
  data("E_LCE")
  SRS <- srs(E_LCE, 0.8, 5)
  expect_equal(sum(!diag(SRS) == 1), 0)
  expect_true(abs(sum(SRS) - 562.5772) < 0.001)
})

test_that("Check srs with wrong inputs",{
  expect_error(srs(E_LCE))
  expect_error(srs(E_LCE, 1.2))
  expect_error(srs(E_LCE, -9))
  expect_error(srs(c(1, 1, 2, 3), 0.8))
  expect_error(srs(E_LCE, dc = 0.7, R = -6))
})

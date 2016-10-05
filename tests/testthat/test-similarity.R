
context("Similarity matrices")

data("E_LCE")

test_that("Check srs works", {
  data("E_LCE")
  SRS <- srs(E_LCE, 0.8, 5)
  expect_equal(sum(!diag(SRS) == 1), 0)
  expect_true(abs(sum(SRS) - 562.5772) < 0.001)
})

test_that("Error in srs with wrong inputs", {
  expect_error(srs(E_LCE))
  expect_error(srs(E_LCE, 1.2))
  expect_error(srs(E_LCE, -9))
  expect_error(srs(c(1, 1, 2, 3), 0.8))
  expect_error(srs(E_LCE, dc = 0.7, R = -6))
})

test_that("Check asrs", {
  ASRS <- asrs(E_LCE, 0.8)
  expect_equal(sum(!diag(ASRS) == 1), 0)
  expect_true(abs(sum(ASRS) - 1456) < 1)
})

test_that("Error in asrs with wrong inputs", {
  expect_error(asrs(E_LCE))
  expect_error(asrs(E_LCE, 1.2))
  expect_error(asrs(E_LCE, -9))
  expect_error(asrs(c(1, 1, 2, 3), 0.8))
})

test_that("Check cts works", {
  CTS <- cts(E_LCE, 0.8)
  expect_equal(sum(!diag(CTS) == 1), 0)
  expect_true(abs(sum(CTS) - 4076.2) < 0.1)
})

test_that("Error in cts with wrong inputs", {
  expect_error(cts(E_LCE))
  expect_error(cts(E_LCE, 1.2))
  expect_error(cts(E_LCE, -9))
  expect_error(cts(c(1, 1, 2, 3), 0.8))
})

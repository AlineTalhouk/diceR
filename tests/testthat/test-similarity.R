context("Similarity matrices")

set.seed(1)
E <- matrix(rep(sample(1:4, 500, replace = TRUE)), nrow = 100)
dc <- 0.8
K <- 4
R <- 5

test_that("Check srs works", {
  SRS <- srs(E = E, dc = dc, R = 5)
  expect_equal(sum(!diag(SRS) == 1), 0)
  expect_equal(abs(sum(SRS)), 813.6889, tolerance = 0.1)
})

test_that("Error in srs with wrong inputs", {
  expect_error(srs(E = E))
  expect_error(srs(E = E, dc = 1.2))
  expect_error(srs(E = E, dc = -9))
  expect_error(srs(E = c(1, 1, 2, 3), dc = 0.8))
  expect_error(srs(E = E, dc = 0.7, R = -6))
})

test_that("Check asrs", {
  ASRS <- asrs(E = E, dc = 0.8)
  expect_equal(sum(!diag(ASRS) == 1), 0)
  expect_equal(abs(sum(ASRS)), 6260.9, tolerance = 0.1)
})

test_that("Error in asrs with wrong inputs", {
  expect_error(asrs(E = E))
  expect_error(asrs(E = E, dc = 1.2))
  expect_error(asrs(E = E, dc = -9))
  expect_error(asrs(E = c(1, 1, 2, 3), dc = 0.8))
})

test_that("Check cts works", {
  CTS <- cts(E = E, dc = 0.8)
  expect_equal(sum(!diag(CTS) == 1), 0)
  expect_equal(abs(sum(CTS)), 7769.9, tolerance = 0.1)
})

test_that("Error in cts with wrong inputs", {
  expect_error(cts(E = E))
  expect_error(cts(E = E, dc = 1.2))
  expect_error(cts(E = E, dc = -9))
  expect_error(cts(E = c(1, 1, 2, 3), dc = 0.8))
})

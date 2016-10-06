
context("External validity indices")

data("FGT", "FGD", "LT", "LD", "E_LCE",
     "final_c1_valid_RandIndex", "final_c2_valid_RandIndex")

test_that("Check ev_accuracy works with basic inputs", {
  expect_equal(ev_accuracy(c(1, 2, 1, 1), c(1, 1, 1, 1)), 1)
  expect_equal(ev_accuracy(c(4, 2, 4, 4), c(4, 4, 4, 4)), 1)
  expect_equal(ev_accuracy(c(4, 4, 4, 4), c(4, 2, 4, 4)), 0.75)
  expect_equal(ev_accuracy(c(4, 4, 4, 4, 8), c(4, 2, 4, 4, 6)), 0.8)
})

test_that("Check ev_accuracy throws error with wrong inputs", {
  expect_error(ev_accuracy(c(1, 2, 3, 4)))
  expect_error(ev_accuracy("a"))
  expect_error(ev_accuracy(c("a", "b", "c", "d")))
  expect_equal(ev_accuracy(c("a", "a", "b"), c("a", "b", "c")), 2 / 3)
})

test_that("Check ev_rand with E_LCE and in case where AR is non-zero", {
  expect_true(abs(ev_rand(E_LCE[, 1], E_LCE[, 2])$AR - 0.2732) <= 0.0001)
  expect_true(abs(ev_rand(E_LCE[, 1], E_LCE[, 2])$HI - 0.3446) <= 0.0001)
  expect_true(abs(ev_rand(E_LCE[, 1], E_LCE[, 2])$MI - 0.3277) <= 0.0001)
  expect_true(abs(ev_rand(E_LCE[, 1], E_LCE[, 2])$RI - 0.6723) <= 0.0001)
  expect_true(ev_rand(final_c1_valid_RandIndex, final_c2_valid_RandIndex)$AR ==
                0)
  expect_true(abs(
    ev_rand(final_c1_valid_RandIndex, final_c2_valid_RandIndex)$RI -
      0.6254
  ) <= 0.0001)
  expect_true(abs(
    ev_rand(final_c1_valid_RandIndex, final_c2_valid_RandIndex)$HI -
      0.2508
  ) <= 0.0001)
  expect_true(abs(
    ev_rand(final_c1_valid_RandIndex, final_c2_valid_RandIndex)$MI -
      0.3746
  ) <= 0.0001)
})

test_that("Check ev_rand with wrong inputs", {
  expect_error(ev_rand(1, 2))
  expect_error(ev_rand(
    data.frame(name = letters[1:4], vals = 1:4),
    data.frame(name = letters[1:4], vals = 1:4)
  ))
  expect_error(ev_rand(c(1, 2, 3, 4),
                       data.frame(name = letters[1:4], vals = 1:4)))
  expect_error(ev_rand(data.frame(name = letters[1:4], vals = 1:4),
                       c(6, 8, 3, 8)))
  expect_error(ev_rand(E_LCE[1:10, 1], E_LCE[, 2]))
  expect_error(ev_rand(matrix(c(1, 2, 3), ncol = 3),
                       matrix(c(1, 2, 3), ncol = 3)))
})

test_that("Check ev_rand with artificial input vectors in case where AR is 0", {
  expect_true(ev_rand(final_c1_valid_RandIndex,
                      final_c2_valid_RandIndex)$AR == 0)
  expect_true(abs(
    ev_rand(final_c1_valid_RandIndex, final_c2_valid_RandIndex)$RI - 0.6254)
    <= 0.0001)
  expect_true(abs(
    ev_rand(final_c1_valid_RandIndex, final_c2_valid_RandIndex)$HI - 0.2508)
    <= 0.0001)
  expect_true(abs(
    ev_rand(final_c1_valid_RandIndex, final_c2_valid_RandIndex)$MI - 0.3746)
    <= 0.0001)
})
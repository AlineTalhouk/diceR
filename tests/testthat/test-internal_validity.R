
context("Internal validity indices")

set.seed(1)
MD <- as.data.frame(matrix(runif(1000, -10, 10), nrow = 100))
set.seed(1)
MT <- sample(1:4, 100, replace = TRUE)

test_that("Check iv_compactness", {
  expect_true(abs(iv_compactness(MD, MT) - 24.1316) < 0.0001)
  MT_k5 <- MT
  MT_k5[10] <- 5
  MT[which(MT == 1)] <- 3
  expect_true(abs(iv_compactness(MD, MT) - 24.5493) < 0.0001)
  expect_true(abs(iv_compactness(MD, MT_k5) - 23.8454) < 0.0001)
})

test_that("iv_compactness throws error with wrong inputs", {
  expect_error(iv_compactness(MD, NULL))
  expect_error(iv_compactness(NULL,c(1,2,3)))
  expect_error(iv_compactness(c(1, 2, 3, 4), c(1, 3, 3, 1)))
})

test_that("PAC can have different bounds", {
  set.seed(1)
  x <- replicate(100, rbinom(100, 4, 0.2))
  y <- consensus_matrix(x)
  expect_error(PAC(y, lower = 0.3, upper = 0.7), NA)
})

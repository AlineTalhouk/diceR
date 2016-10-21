
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

test_that("Check iv_db_dunn with MD and MT", {
  expect_true(abs(iv_db_dunn(MD, MT)$DB - 4.6724) <= 0.001)
  expect_true(abs(iv_db_dunn(MD, MT)$Dunn - 0.3734) <= 0.001)
  expect_identical(iv_db_dunn(MD, MT), iv_db_dunn(MD, data.frame(MT)))
})

test_that("Check iv_db_dunn with wrong inputs", {
  expect_error(iv_db_dunn(MD, MT[1:99]))
  expect_error(iv_db_dunn(c(1, 2, 3, 4), c(1, 2, 3, 4)))
  expect_error(iv_db_dunn(MD, matrix(c(1, 2, 3), ncol = 3)))
})

test_that("PAC can have different bounds", {
  set.seed(1)
  x <- replicate(100, rbinom(100, 4, 0.2))
  y <- consensus_matrix(x)
  expect_error(PAC(y, lower = 0.3, upper = 0.7), NA)
})

test_that("CHI wrapper same as original", {
  expect_identical(iv_chi(MD, MT), clusterSim::index.G1(MD, MT))
})
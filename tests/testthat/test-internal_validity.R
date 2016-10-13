
context("Internal validity indices")

data("FGT", "FGD", "LT", "LD", "E_LCE")

test_that("Check iv_compactness", {
  labels <- E_LCE[, 1]
  expect_true(abs(iv_compactness(FGD, labels) - 5.8559) < 0.0001)
  labels[which(labels == 1)] <- 3
  labels[1] <- 1
  expect_true(abs(iv_compactness(FGD, labels) - 6.2886) < 0.0001)
})

test_that("iv_compactness throws error with wrong inputs", {
  expect_error(iv_compactness(FGD, E_LCE[-18, 1]))
  expect_error(iv_compactness(c(1, 2, 3, 4), c(1, 3, 3, 1)))
})

test_that("Check iv_db_dunn with FGD and FGT", {
  expect_true(abs(iv_db_dunn(FGD, FGT)$DB - 1.1441) <= 0.001)
  expect_true(abs(iv_db_dunn(FGD, FGT)$Dunn - 1.3570) <= 0.001)
  expect_true(abs(iv_db_dunn(LD, LT)$DB - 3.2857) <= 0.001)
  expect_true(abs(iv_db_dunn(LD, LT)$Dunn - 0.5939) <= 0.001)
})

test_that("Check iv_db_dunn with wrong inputs", {
  expect_error(iv_db_dunn(FGD, FGT[1:99]))
  expect_error(iv_db_dunn(c(1, 2, 3, 4), c(1, 2, 3, 4)))
  expect_error(iv_db_dunn(FGD, matrix(c(1, 2, 3), ncol = 3)))
})

test_that("Check iv_sumsq with FGD and FGT", {
  test_FGD <- iv_sumsq(data = FGD, labels = FGT, k = 4)
  expect_true(nrow(test_FGD$Tot) == 12)
  expect_true(ncol(test_FGD$Tot) == 12)
  expect_true(abs(sum(test_FGD$Tot) - 2622.8) <= 0.1)
  expect_true(abs(sum(test_FGD$Sintra) - 15.0846) <= 0.1)
  expect_true(nrow(test_FGD$Sinter) == 4)
  expect_true(ncol(test_FGD$Sinter) == 4)
  expect_true(nrow(test_FGD$W) == 12)
  expect_true(ncol(test_FGD$W) == 12)
  expect_true(abs(sum(test_FGD$Sinter) - 80.3340) <= 0.0001)
  expect_true(abs(sum(test_FGD$W) - 1528.8) <= 0.1)

})

test_that("Check iv_sumsq with LD and LT with k=4", {
  test_LD <- iv_sumsq(data = LD, labels = LT, k = 4)
  expect_true(nrow(test_LD$Tot) == 1081)
  expect_true(ncol(test_LD$Tot) == 1081)
  expect_true(abs(sum(test_LD$Tot) - 1653379.08) <= 0.01)
  expect_true(abs(sum(test_LD$W) - 1613148.64) <= 0.01)
  expect_true(nrow(test_LD$W) == 1081)
  expect_true(ncol(test_LD$W) == 1081)
  expect_true(nrow(test_LD$Sinter) == 4)
  expect_true(ncol(test_LD$Sinter) == 4)
  expect_true(abs(sum(test_LD$Sintra) - 61.17955) <= 0.00001)
  expect_true(abs(sum(test_LD$Sinter) - 111.7192) <= 0.0001)
  expect_true(nrow(test_LD$B) == 1081)
  expect_true(ncol(test_LD$B) == 1081)
  expect_true(abs(sum(test_LD$B) - 40230.439) <= 0.001)
})

test_that("Check iv_sumsq with LD and LT with k=1", {
  test_LD_k1 <- iv_sumsq(LD, LT, 1)
  expect_true(sum(test_LD_k1$B) == 0)
  expect_true(sum(test_LD_k1$Sinter) == 0)
  expect_true(abs(sum(test_LD_k1$Sintra) - 31.9592) <= 0.0001)
  expect_true(abs(sum(test_LD_k1$W) - 1653379.08) <= 0.01)
  expect_true(abs(sum(test_LD_k1$Tot) - 1653379.08) <= 0.01)
})

test_that("Check iv_sumsq with wrong inputs", {
  expect_error(iv_sumsq(LD, FGT, 4))
  expect_error(iv_sumsq(FGD, FGT, "a"))
  expect_error(iv_sumsq(FGD, sample(
    1:4, size = 18, replace = TRUE
  ), 4))
  expect_error(iv_sumsq(FGD, as.factor(FGT), 4))
  expect_error(iv_sumsq(
    data = as.factor(FGD[, 1]),
    labels = FGT,
    k = 4
  ))
})

test_that("PAC can have different bounds", {
  set.seed(1)
  x <- replicate(100, rbinom(100, 4, 0.2))
  y <- consensus_matrix(x)
  expect_error(PAC(y, lower = 0.3, upper = 0.7), NA)
})

context("validility scores")

data("FGT")
data("FGD")
data("LT")
data("LD")
data("E_LCE")

test_that("Check valid_CA works with basic inputs", {
  expect_equal(valid_CA(c(1, 2, 1, 1), c(1, 1, 1, 1)), 1)
  expect_equal(valid_CA(c(4, 2, 4, 4), c(4, 4, 4, 4)), 1)
  expect_equal(valid_CA(c(4, 4, 4, 4), c(4, 2, 4, 4)), 0.75)
  expect_equal(valid_CA(c(4, 4, 4, 4, 8), c(4, 2, 4, 4, 6)), 0.8)
})

test_that("Check valid_CA throws error with wrong inputs", {
  expect_error(valid_CA(c(1, 2, 3, 4)))
  expect_error(valid_CA("a"))
  expect_error(valid_CA(c("a", "b", "c", "d")))
  expect_equal(valid_CA(c("a", "a", "b"), c("a", "b", "c")), 2 / 3)
})

test_that("Check valid_compactness", {
  labels <- E_LCE[, 1]
  expect_true(abs(valid_compactness(FGD, labels) - 5.8559) < 0.0001)
  labels[which(labels == 1)] <- 3
  labels[1] <- 1
  expect_true(abs(valid_compactness(FGD, labels) - 6.2886) < 0.0001)
})

test_that("valid_compactness throws error with wrong inputs", {
  expect_error(valid_compactness(FGD, E_LCE[-18, 1]))
  expect_error(valid_compactness(c(1, 2, 3, 4), c(1, 3, 3, 1)))
})


test_that("Check valid_DbDunn with FGD and FGT", {
  expect_true(abs(valid_DbDunn(FGD, FGT)$DB - 1.1441) <= 0.001)
  expect_true(abs(valid_DbDunn(FGD, FGT)$Dunn - 1.3570) <= 0.001)
  expect_true(abs(valid_DbDunn(LD, LT)$DB - 3.2857) <= 0.001)
  expect_true(abs(valid_DbDunn(LD, LT)$Dunn - 0.5939) <= 0.001)

})

test_that("Check valid_DbDunn with wrong inputs", {
  expect_error(valid_DbDunn(FGD, FGT[1:99]))
  expect_error(valid_DbDunn(c(1, 2, 3, 4), c(1, 2, 3, 4)))
  expect_error(valid_DbDunn(FGD,matrix(c(1,2,3),ncol=3)))
})


data("E_LCE")
data("final_c1_valid_RandIndex")
data("final_c2_valid_RandIndex")

test_that("Check valid_RandIndex with E_LCE and in case where AR is non-zero", {
  expect_true(abs(valid_RandIndex(E_LCE[, 1], E_LCE[, 2])$AR - 0.2732) <= 0.0001)
  expect_true(abs(valid_RandIndex(E_LCE[, 1], E_LCE[, 2])$HI - 0.3446) <= 0.0001)
  expect_true(abs(valid_RandIndex(E_LCE[, 1], E_LCE[, 2])$MI - 0.3277) <= 0.0001)
  expect_true(abs(valid_RandIndex(E_LCE[, 1], E_LCE[, 2])$RI - 0.6723) <= 0.0001)
  expect_true(valid_RandIndex(final_c1_valid_RandIndex, final_c2_valid_RandIndex)$AR ==
                0)
  expect_true(abs(
    valid_RandIndex(final_c1_valid_RandIndex, final_c2_valid_RandIndex)$RI -
      0.6254
  ) <= 0.0001)
  expect_true(abs(
    valid_RandIndex(final_c1_valid_RandIndex, final_c2_valid_RandIndex)$HI -
      0.2508
  ) <= 0.0001)
  expect_true(abs(
    valid_RandIndex(final_c1_valid_RandIndex, final_c2_valid_RandIndex)$MI -
      0.3746
  ) <= 0.0001)
})

test_that("Check valid_RandIndex with wrong inputs", {
  expect_error(valid_RandIndex(1, 2))
  expect_error(valid_RandIndex(
    data.frame(name = letters[1:4], vals = 1:4),
    data.frame(name = letters[1:4], vals = 1:4)
  ))
  expect_error(valid_RandIndex(c(1, 2, 3, 4), data.frame(name = letters[1:4], vals =
                                                           1:4)))
  expect_error(valid_RandIndex(data.frame(name = letters[1:4], vals = 1:4), c(6, 8, 3, 8)))
  expect_error(valid_RandIndex(E_LCE[1:10, 1], E_LCE[, 2]))
  expect_error(valid_RandIndex(matrix(c(1,2,3),ncol=3),matrix(c(1,2,3),ncol=3)))
})

test_that("Check valid_RandIndex with artificial input vectors in case where AR is 0",
          {
            expect_true(valid_RandIndex(final_c1_valid_RandIndex, final_c2_valid_RandIndex)$AR ==
                          0)
            expect_true(abs(
              valid_RandIndex(final_c1_valid_RandIndex, final_c2_valid_RandIndex)$RI -
                0.6254
            ) <= 0.0001)
            expect_true(abs(
              valid_RandIndex(final_c1_valid_RandIndex, final_c2_valid_RandIndex)$HI -
                0.2508
            ) <= 0.0001)
            expect_true(abs(
              valid_RandIndex(final_c1_valid_RandIndex, final_c2_valid_RandIndex)$MI -
                0.3746
            ) <= 0.0001)
          })

data("LD")
data("LT")
data("FGD")
data("FGT")

test_that("Check valid_sumsqures with FGD and FGT", {
  #Test using FGD and FGT
  test_FGD <- valid_sumsqures(data = FGD, labels = FGT, k = 4)
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

test_that("Check valid_sumsqures with LD and LT with k=4", {
  test_LD <- valid_sumsqures(data = LD, labels = LT, k = 4)
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

test_that("Check valid_sumsqures with LD and LT with k=1", {
  #Test using LD and LT with k=1
  test_LD_k1 <- valid_sumsqures(LD, LT, 1)
  expect_true(sum(test_LD_k1$B) == 0)
  expect_true(sum(test_LD_k1$Sinter) == 0)
  expect_true(abs(sum(test_LD_k1$Sintra) - 31.9592) <= 0.0001)
  expect_true(abs(sum(test_LD_k1$W) - 1653379.08) <= 0.01)
  expect_true(abs(sum(test_LD_k1$Tot) - 1653379.08) <= 0.01)
})

test_that("Check valid_sumsqures with wrong inputs", {
  #Test with error cases
  expect_error(valid_sumsqures(LD, FGT, 4))
  expect_error(valid_sumsqures(FGD, FGT, "a"))
  expect_error(valid_sumsqures(FGD, sample(
    1:4, size = 18, replace = TRUE
  ), 4))
  expect_error(valid_sumsqures(FGD, as.factor(FGT), 4))
  expect_error(valid_sumsqures(
    data = as.factor(FGD[, 1]),
    labels = FGT,
    k = 4
  ))
})

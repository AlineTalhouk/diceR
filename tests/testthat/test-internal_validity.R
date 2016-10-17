
context("Internal validity indices")

set.seed(1)
MD<-as.data.frame(matrix(runif(1000,-10,10),nrow=100,byrow=FALSE))
set.seed(1)
MT<-sample(1:4,100,replace=TRUE)

test_that("Check iv_compactness", {
  expect_true(abs(iv_compactness(MD, MT) - 24.1316) < 0.0001)
  MT_k5<-MT
  MT_k5[10]<-5
  MT[which(MT == 1)] <- 3
  expect_true(abs(iv_compactness(MD, MT) - 24.5493) < 0.0001)
  expect_true(abs(iv_compactness(MD,MT_k5)-23.8454)<0.0001)
})

test_that("iv_compactness throws error with wrong inputs", {
  expect_error(iv_compactness(MD, NULL))
  expect_error(iv_compactness(NULL,c(1,2,3)))
  expect_error(iv_compactness(c(1, 2, 3, 4), c(1, 3, 3, 1)))
})

test_that("Check iv_db_dunn with MD and MT", {
  set.seed(1)
  MD<-as.data.frame(matrix(runif(1000,-10,10),nrow=100,byrow=FALSE))
  set.seed(1)
  MT<-sample(1:4,100,replace=TRUE)
  expect_true(abs(iv_db_dunn(MD, MT)$DB - 4.6724) <= 0.001)
  expect_true(abs(iv_db_dunn(MD, MT)$Dunn - 0.3734) <= 0.001)
})

test_that("Check iv_db_dunn with wrong inputs", {
  expect_error(iv_db_dunn(MD, MT[1:99]))
  expect_error(iv_db_dunn(c(1, 2, 3, 4), c(1, 2, 3, 4)))
  expect_error(iv_db_dunn(MD, matrix(c(1, 2, 3), ncol = 3)))
})

test_that("Check iv_sumsq with MD and MT", {
  set.seed(1)
  MD<-as.data.frame(matrix(runif(1000,-10,10),nrow=100,byrow=FALSE))
  set.seed(1)
  MT<-sample(1:4,100,replace=TRUE)
  test_MD <- iv_sumsq(data = MD, labels = MT, k = 4)
  expect_true(nrow(test_MD$Tot) == 10)
  expect_true(ncol(test_MD$Tot) == 10)
  expect_true(abs(sum(test_MD$Tot) - 40863) < 1)
  expect_true(abs(sum(test_MD$Sintra) - 67.6140) < 0.0001)
  expect_true(nrow(test_MD$Sinter) == 4)
  expect_true(ncol(test_MD$Sinter) == 4)
  expect_true(nrow(test_MD$W) == 10)
  expect_true(ncol(test_MD$W) == 10)
  expect_true(abs(sum(test_MD$Sinter) - 118.5502) < 0.0001)
  expect_true(abs(sum(test_MD$W) - 35634) < 1)

})


test_that("Check iv_sumsq with MD and MT with k=1", {
  set.seed(1)
  MD<-as.data.frame(matrix(runif(1000,-10,10),nrow=100,byrow=FALSE))
  set.seed(1)
  MT<-sample(1:4,100,replace=TRUE)
  test_MD_k1 <- iv_sumsq(MD, MT, 1)
  expect_true(sum(test_MD_k1$B) == 0)
  expect_true(sum(test_MD_k1$Sinter) == 0)
  expect_true(abs(sum(test_MD_k1$Sintra) - 17.9818) < 0.0001)
  expect_true(abs(sum(test_MD_k1$W) - 40863) < 1)
  expect_true(abs(sum(test_MD_k1$Tot) - 40863) < 1)
})

test_that("Check iv_sumsq with wrong inputs", {
  expect_error(iv_sumsq(NULL, MT, 4))
  expect_error(iv_sumsq(MD, MT, "a"))
  expect_error(iv_sumsq(MD, sample(
    1:4, size = 18, replace = TRUE
  ), 4))
  expect_error(iv_sumsq(MD, as.factor(MT), 4))
  expect_error(iv_sumsq(
    data = as.factor(MD[, 1]),
    labels = MT,
    k = 4
  ))
})

test_that("PAC can have different bounds", {
  set.seed(1)
  x <- replicate(100, rbinom(100, 4, 0.2))
  y <- consensus_matrix(x)
  expect_error(PAC(y, lower = 0.3, upper = 0.7), NA)
})
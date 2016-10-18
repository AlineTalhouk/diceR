
context("ConClustEnsemble")

data(hgsc)
dat <- hgsc[, -1]
rownames(dat) <- hgsc[, 1]
k <- 4

test_that("Check ConClustEnsemble with hcAEucl and kmEucl with hgsc data", {
  reps <- 5
  method <- c("hcAEucl", "kmEucl")
  E <- ConClustEnsemble(X = t(dat), k = k, reps = reps, method = method)
  expect_equal(sum(!complete.cases(E)), 0)
  expect_equal(ncol(dat), nrow(E))
  expect_equal(ncol(E), reps * length(method))
})

test_that("Check ConClustEnsemble with just hcAEucl", {
  reps <- 6
  E <- ConClustEnsemble(X = t(dat), k = k, reps = reps, method = "hcAEucl")
  expect_equal(sum(!complete.cases(E)), 0)
  expect_equal(ncol(dat), nrow(E))
  expect_equal(ncol(E), reps)
})

test_that("Check ConClustEnsemble with all the methods", {
  reps <- 5
  x <- dat[seq_len(200), seq_len(100)]
  E <- ConClustEnsemble(X = t(x), k = k, reps = reps, method = NULL)
  expect_equal(sum(!complete.cases(E)), 0)
  expect_equal(ncol(x), nrow(E))
  expect_equal(ncol(E), reps * 12)
})

test_that("Check ConClustEnsemble with wrong inputs", {
  expect_error(ConClustEnsemble(X = t(dat), k = k, reps = 1,
                                method = c("hcAEucl", "kmEucl")))
  expect_error(ConClustEnsemble(X = c(1, 2, 3), k = k, reps = 5,
                                method = c("hcAEucl", "kmEucl")))
  expect_error(ConClustEnsemble(X = t(dat), k = k, reps = 5,
                                method = c("a", "b", "c")))
})
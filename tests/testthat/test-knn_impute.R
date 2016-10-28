
context("KNN imputation")

test_that("fewer or equal number of NA after imputation", {
  set.seed(2)
  x <- replicate(100, rnorm(100))
  CC <- ConClust(x, nc = 4, reps = 10, method = "scRbf", progress = FALSE)
  orig.cl <- CC[, 1, 1, 1]
  impute.cl <- knn_impute(orig.cl, x)
  expect_lte(sum(is.na(impute.cl)), sum(is.na(orig.cl)))
})
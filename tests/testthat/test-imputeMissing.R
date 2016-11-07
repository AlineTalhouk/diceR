
context("Impute missing values")

data(hgsc)
data <- t(hgsc[, -1])
E <- ConClust(data, nc = 4, reps = 10,
              method = c("hcAEucl", "kmEucl", "scRbf"), progress = FALSE)

test_that("knn has less than or equal number of NA after imputation", {
  E_imputed <- imputeMissing(E, data)
  E_knn <- E_imputed$E_imputed
  E_complete <- E_imputed$E_imputed2
  expect_lte(sum(is.na(E_knn)), sum(is.na(E)))
  expect_equal(sum(is.na(E_complete)), 0)
})

test_that("Check imputeMissing with hgsc data with imputeALL TRUE", {
  E_imputed <- imputeMissing(E, data)
  expect_true(sum(!complete.cases(E_imputed$E_imputed2)) == 0)
  expect_equal(nrow(E_imputed$E_imputed2), dim(E)[1])
  expect_equal(ncol(E_imputed$E_imputed2), dim(E)[2] * dim(E)[3])
})

test_that("Check imputeMissing with hgsc data with imputeALL FALSE", {
  E_imputed <- imputeMissing(E, data, imputeALL = FALSE)
  expect_true(sum(!complete.cases(E_imputed)) < sum(!complete.cases(E)))
  expect_equal(nrow(E_imputed$E_imputed), dim(E)[1])
  expect_equal(ncol(E_imputed$E_imputed), dim(E)[2])
})

test_that("Check imputeMissing throws error with wrong inputs", {
  expect_error(imputeMissing(E, data[-1, ]))
  expect_error(imputeMissing(E = sample(1:489), data = data))
  expect_error(imputeMissing(E = E, data = as.data.frame(data)))
  expect_error(imputeMissing(E = E, data = data, imputeALL = 2))
})
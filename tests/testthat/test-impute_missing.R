
context("Impute missing values")

data(hgsc)
data <- t(hgsc[, -1])
E <- ConClust(data, nk = 4, reps = 10,
              method = c("hcAEucl", "kmEucl", "scRbf"), progress = FALSE)
E_imputed <- impute_missing(E, data)

test_that("knn has less than or equal number of NA after imputation", {
  E_knn <- E_imputed$knn
  expect_lte(sum(is.na(E_knn)), sum(is.na(E)))
})

test_that("adding majority vote elimnates all missing values", {
  E_complete <- E_imputed$complete
  expect_equal(sum(is.na(E_complete)), 0)
  expect_equal(nrow(E_complete), dim(E)[1])
  expect_equal(ncol(E_complete), prod(dim(E)[2], dim(E)[3]))
})

test_that("errors thrown with incorrect inputs", {
  expect_error(impute_missing(E, data[-1, ]))
  expect_error(impute_missing(E = sample(1:489), data = data))
  expect_error(impute_missing(E = E, data = as.data.frame(data)))
})
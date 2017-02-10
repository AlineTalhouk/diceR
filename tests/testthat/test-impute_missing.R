
context("Impute missing values")

data(hgsc)
data <- t(hgsc[, -1])
E <- consensus_cluster(data, nk = 4, reps = 5,
                       algorithms = c("hc", "km", "sc"),
                       progress = FALSE)
E_imputed <- impute_missing(E, data, nk = 4)

test_that("adding majority vote elimnates all missing values", {
  E_complete <- E_imputed
  expect_equal(sum(is.na(E_complete)), 0)
  expect_equal(nrow(E_complete), dim(E)[1])
  expect_equal(ncol(E_complete), prod(dim(E)[2], dim(E)[3]))
})

test_that("errors thrown with incorrect inputs", {
  expect_error(impute_missing(E, data[-1, ]))
  expect_error(impute_missing(E = sample(1:489), data = data))
  expect_error(impute_missing(E = E, data = as.data.frame(data)))
})

test_that("majority completion not performed for a single assignment", {
  E1 <- consensus_cluster(data, nk = 4, reps = 1, algorithms = "hc",
                          progress = FALSE)
  expect_equal(ncol(impute_missing(E1, data, nk = 4)), 1)
})
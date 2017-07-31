context("Impute missing values")

data(hgsc)
hgsc <- hgsc[1:40, 1:30]
E <- consensus_cluster(hgsc, nk = 4, reps = 5,
                       algorithms = c("hc", "km", "sc"), progress = FALSE)
E_imputed <- impute_missing(E, hgsc, nk = 4)

test_that("adding majority vote elimnates all missing values", {
  E_complete <- E_imputed
  expect_equal(sum(is.na(E_complete)), 0)
  expect_equal(nrow(E_complete), dim(E)[1])
  expect_equal(ncol(E_complete), prod(dim(E)[2], dim(E)[3]))
})

test_that("errors thrown with incorrect inputs", {
  expect_error(impute_missing(E, hgsc[-1, ]))
  expect_error(impute_missing(E = sample(1:489), data = hgsc))
  expect_error(impute_missing(E = E, data = as.data.frame(hgsc)))
})

test_that("majority completion not performed for a single assignment", {
  E1 <- consensus_cluster(hgsc, nk = 4, reps = 1, algorithms = "hc",
                          progress = FALSE)
  expect_equal(ncol(impute_missing(E1, hgsc, nk = 4)), 1)
})

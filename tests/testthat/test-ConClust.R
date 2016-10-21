
context("Consensus clustering")

data(hgsc)
dat <- t(hgsc[, -1])

test_that("Output is an array", {
  x <- ConClust(dat, k = 4, reps = 10, method = "hcAEucl", progress = FALSE)
  expect_is(x, "array")
})

test_that("No method given means all methods", {
  expect_identical(ConClust(dat, k = 4, reps = 1, progress = FALSE),
                   ConClust(dat, k = 4, reps = 1, progress = FALSE,
                            parallel = TRUE))
})

test_that("Output can be saved with or without time in file name", {
  x1 <- ConClust(dat, k = 4, reps = 10, method = "hcAEucl", progress = FALSE,
                 save = TRUE)
  x2 <- ConClust(dat, k = 4, reps = 10, method = "hcAEucl", progress = FALSE,
                 save = TRUE, time.saved = TRUE)
  expect_identical(x1, x2)
  file.remove(list.files(pattern = "ConClustOutput"))
})

test_that("Parallel computing can be used", {
  x1 <- ConClust(dat, k = 4, reps = 10, method = "hcAEucl", progress = TRUE, parallel = TRUE)
  x2 <- ConClust(dat, k = 4, reps = 10, method = "hcAEucl", progress = FALSE, parallel = TRUE)
  x3 <- ConClust(dat, k = 4, reps = 10, method = "hcAEucl", progress = TRUE, parallel = FALSE)
  expect_identical(x1, x2)
  expect_identical(x2, x3)
})

test_that("Parallel computing automatically used under certain settings", {
  x3 <- ConClust(dat, k = 4, reps = 100, method = "hcAEucl", progress = TRUE)
  expect_error(x3, NA)
})

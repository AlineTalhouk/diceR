
context("Consensus clustering")

data(hgsc)
dat <- t(hgsc[, -1])

test_that("Output is an array", {
  x <- ConClust(dat, nc = 2:4, reps = 10, method = "hcAEucl", progress = FALSE)
  expect_is(x, "array")
})

test_that("No method given means all methods", {
  x1 <- ConClust(dat, nc = 2:4, reps = 1, progress = FALSE)
  x2 <- ConClust(dat, nc = 2:4, reps = 1, progress = FALSE, parallel = TRUE)
  closeAllConnections()
  expect_identical(x1, x2)
})

test_that("Output can be saved with or without time in file name", {
  x1 <- ConClust(dat, nc = 2:4, reps = 10, method = "hcAEucl", progress = FALSE,
                 save = TRUE)
  x2 <- ConClust(dat, nc = 2:4, reps = 10, method = "hcAEucl", progress = FALSE,
                 save = TRUE, time.saved = TRUE)
  expect_identical(x1, x2)
  file.remove(list.files(pattern = "ConClustOutput"))
})

test_that("Parallel computing can be used", {
  x1 <- ConClust(dat, nc = 2:4, reps = 10, method = "hcAEucl", progress = TRUE,
                 parallel = TRUE, ncores = 4)
  x2 <- ConClust(dat, nc = 2:4, reps = 10, method = "hcAEucl", progress = FALSE,
                 parallel = TRUE, ncores = 4)
  x3 <- ConClust(dat, nc = 2:4, reps = 10, method = "hcAEucl", progress = TRUE,
                 parallel = FALSE)
  closeAllConnections()
  expect_identical(x1, x2)
  expect_identical(x2, x3)
})

test_that("Parallel computing automatically used under certain settings", {
  expect_error(ConClust(dat, nc = 2:4, reps = 100, method = "hcAEucl", progress = TRUE), NA)
})

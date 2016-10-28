
context("Consensus clustering")

data(hgsc)
dat <- t(hgsc[, -1])

test_that("Output is an array", {
  x <- ConClust(dat, nc = 2:4, reps = 10, method = "hcAEucl", progress = FALSE)
  expect_is(x, "array")
})

test_that("No method given means all methods", {
  x1 <- ConClust(dat, nc = 2:4, reps = 1, progress = FALSE)
  expect_error(x1, NA)
})

test_that("Output can be saved with or without time in file name", {
  x1 <- ConClust(dat, nc = 2:4, reps = 10, method = "hcAEucl", progress = FALSE,
                 save = TRUE)
  x2 <- ConClust(dat, nc = 2:4, reps = 10, method = "hcAEucl", progress = FALSE,
                 save = TRUE, time.saved = TRUE)
  expect_identical(x1, x2)
  file.remove(list.files(pattern = "ConClustOutput"))
})
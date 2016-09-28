
context("Consensus clustering")

data(hgsc)

test_that("Output is an array", {
  x <- ConClust(hgsc, k = 4, reps = 10, method = "hcAEucl", save = FALSE)
  expect_is(x, "array")
})

test_that("No method given means all methods", {
  x <- ConClust(hgsc, k = 4, reps = 1, save = FALSE)
  expect_equal(dim(x)[3], 12)
})

test_that("Output can be saved with or without time in file name", {
  x1 <- ConClust(hgsc, k = 4, reps = 10, method = "hcAEucl", save = TRUE)
  x2 <- ConClust(hgsc, k = 4, reps = 10, method = "hcAEucl", save = TRUE,
                 time.saved = TRUE)
  expect_identical(x1, x2)
  file.remove(list.files(pattern = "ConClustOutput"))
})

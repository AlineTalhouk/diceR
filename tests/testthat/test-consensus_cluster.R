
context("Consensus clustering")

data(hgsc)
dat <- t(hgsc[, -1])

test_that("Output is an array", {
  x <- consensus_cluster(dat, nk = 2:4, reps = 5, algorithms = "hc",
                         progress = FALSE)
  expect_is(x, "array")
})

test_that("No algorithms given means all algorithmss", {
  x1 <- consensus_cluster(dat, nk = 4, reps = 1, progress = FALSE)
  expect_error(x1, NA)
})

test_that("Output can be saved with or without time in file name", {
  x1 <- consensus_cluster(dat, nk = 2:4, reps = 5, algorithms = "hc",
                          progress = FALSE, save = TRUE)
  x2 <- consensus_cluster(dat, nk = 2:4, reps = 5, algorithms = "hc",
                          progress = FALSE, save = TRUE, time.saved = TRUE)
  expect_identical(x1, x2)
  file.remove(list.files(pattern = "CCOutput"))
})
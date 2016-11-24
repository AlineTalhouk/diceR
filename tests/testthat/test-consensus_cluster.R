
context("Consensus clustering")

data(hgsc)
dat <- t(hgsc[, -1])

test_that("No algorithms means all algorithms, output is an array", {
  x1 <- consensus_cluster(dat, nk = 4, reps = 1, progress = FALSE)
  expect_error(x1, NA)
  expect_is(x1, "array")
})

test_that("Output can be saved with or without time in file name", {
  x1 <- consensus_cluster(dat, nk = 2:4, reps = 5, algorithms = "hc",
                          progress = FALSE, save = TRUE)
  x2 <- consensus_cluster(dat, nk = 2:4, reps = 5, algorithms = "hc",
                          progress = FALSE, save = TRUE, time.saved = TRUE)
  expect_identical(x1, x2)
  file.remove(list.files(pattern = "CCOutput"))
})

test_that("Progress bar increments across entire function call", {
  assign("my_dist", function(x) stats::dist(x, method = "manhattan"), pos = 1)
  x3 <- consensus_cluster(dat, nk = 2, reps = 5, algorithms = c("nmf", "hc", "ap"),
                          distance = c("spear", "my_dist") , nmf.method = "lee",
                          progress = TRUE)
  expect_error(x3, NA)
})

test_that("Able to call only spearman distance", {
  x4 <- consensus_cluster(dat, nk = 2, reps = 5, algorithms = "hc",
                          distance = "spear")
  expect_error(x4, NA)
})
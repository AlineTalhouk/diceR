
context("Consensus summary")

data(hgsc)
dat <- t(hgsc[, -1])
x <- ConClust(dat, k = 4, reps = 10, method = "hcAEucl", progress = FALSE)
y1 <- consensus_summary(x, k = 2, progress = FALSE)
y2 <- consensus_summary(x, k = 2, progress = FALSE, save = TRUE)

test_that("Summary is a list of two elements", {
  expect_is(y1, "list")
  expect_length(y1$hcAEucl, 2)
})

test_that("Summary can be saved", {
  expect_identical(y1, y2)
  file.remove(list.files(pattern = "results_CC"))
})

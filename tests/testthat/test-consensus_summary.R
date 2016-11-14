
context("Consensus summary")

data(hgsc)
dat <- t(hgsc[, -1])
x <- ConClust(dat, nk = 2:4, reps = 10, method = "hcAEucl", progress = FALSE)
y <- consensus_summary(x)

test_that("Summary is a list of layered elements", {
  expect_is(y, "list")
  expect_length(y$`2`$hcAEucl, 2)
})
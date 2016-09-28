
context("Consensus summary")

data(hgsc)
x <- ConClust(hgsc, k = 4, reps = 10, method = "hcAEucl", save = FALSE)
y1 <- consensus_summary(x, k = 2)
y2 <- consensus_summary(x, k = 2, save = TRUE)

test_that("Summary is a list of two elements", {
  expect_is(y1, "list")
  expect_length(y1$hcAEucl, 2)
})

test_that("Summary can be saved", {
  expect_identical(y1, y2)
  file.remove(list.files(pattern = "results_CC"))
})

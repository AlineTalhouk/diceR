
context("Consensus combine, evaluate, and weigh")

set.seed(911)
x <- matrix(rnorm(1000), nrow = 10)
CC1 <- ConClust(x, k = 4, reps = 5, method = "apEucl", progress = FALSE)
CC2 <- ConClust(x, k = 4, reps = 5, method = "gmmBIC", progress = FALSE)
CC1.summ <- consensus_summary(CC1, k = 4, progress = FALSE)
CC2.summ <- consensus_summary(CC2, k = 4, progress = FALSE)

an <- c("A1", "A2")
y1 <- consensus_combine(CC1.summ, CC2.summ, element = "matrix")
y2 <- consensus_combine(CC1.summ, CC2.summ, element = "class")
y3 <- consensus_combine(CC1.summ, CC2.summ, element = "matrix", alg.names = an)
y4 <- consensus_combine(CC1.summ, CC2.summ, element = "class", alg.names = an)

z1 <- consensus_evaluate(x, cl.mat = y2, cons.mat = y1)
z2 <- consensus_evaluate(x, cl.mat = y2, cons.mat = y1, alg.names =  an)

test_that("combining results has expected lengths", {
  expect_length(y1, dim(CC1)[3] + dim(CC2)[3])
  expect_equal(ncol(y2), dim(CC1)[3] + dim(CC2)[3])
})

test_that("names can be overwritten", {
  expect_identical(names(y3), an)
  expect_identical(colnames(y4), an)
})

test_that("comparing results works", {
  expect_error(z1, NA)
  expect_error(z2, NA)
})

test_that("weighing works", {
  expect_error(consensus_weigh(z1$internal), NA)
})

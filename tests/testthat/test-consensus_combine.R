
context("Consensus combine, evaluate, and weigh")

set.seed(911)
x <- matrix(rnorm(1000), nrow = 10)
CC1 <- ConClust(x, nc = 2:4, reps = 5, method = "apEucl", progress = FALSE)
CC2 <- ConClust(x, nc = 2:4, reps = 5, method = "gmmBIC", progress = FALSE)
CC1.summ <- consensus_summary(CC1, progress = FALSE)
CC2.summ <- consensus_summary(CC2, progress = FALSE)

an <- c("A1", "A2")
y1 <- consensus_combine(CC1.summ, CC2.summ, element = "matrix")
y2 <- consensus_combine(CC1.summ, CC2.summ, element = "class")
y3 <- consensus_combine(CC1.summ, CC2.summ, element = "matrix", alg.names = an)
y4 <- consensus_combine(CC1.summ, CC2.summ, element = "class", alg.names = an)

z1a <- consensus_evaluate(x, cl.mat = y2, cons.mat = y1)
z1b <- consensus_evaluate(x, cl.mat = y2, cons.mat = y1, alg.names =  an)
set.seed(1)
ref.cl <- sample(1:4, 10, replace = TRUE)
z2 <- consensus_evaluate(x, cl.mat = y2, cons.mat = y1, ref.cl)

test_that("combining results has expected lengths", {
  expect_length(y1, dim(CC1)[3] + dim(CC2)[3])
  expect_equal(ncol(y2), dim(CC1)[3] + dim(CC2)[3])
})

test_that("names can be overwritten", {
  expect_identical(names(y3), an)
  expect_identical(colnames(y4), an)
})

test_that("comparing results works", {
  expect_error(z1a, NA)
  expect_error(z1b, NA)
  expect_length(z1a, 1)
  expect_length(z2, 2)
})

test_that("weighing works", {
  expect_error(consensus_weigh(z1a$internal), NA)
})

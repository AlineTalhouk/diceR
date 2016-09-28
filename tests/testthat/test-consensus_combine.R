
context("Consensus results combine and compare")

set.seed(911)
x <- matrix(rnorm(1000), nrow = 10)
CC1 <- ConClust(x, k = 4, reps = 5, method = "apEucl", save = FALSE)
CC2 <- ConClust(x, k = 4, reps = 5, method = "gmmBIC", save = FALSE)
CC1.summ <- consensus_summary(CC1, k = 4)
CC2.summ <- consensus_summary(CC2, k = 4)
CCP <- ConClustPlus(x, k = 4, reps = 5, save = FALSE)

y1 <- consensus_combine(CC1.summ, CC2.summ, res.CCP = CCP, k = 4,
                        element = "matrix")
y2 <- consensus_combine(CC1.summ, CC2.summ, res.CCP = CCP, k = 4,
                        element = "class")
y3 <- consensus_combine(CC1.summ, CC2.summ, res.CCP = CCP, k = 4,
                        element = "matrix", alg.names = paste0("A", 1:8))
y4 <- consensus_combine(CC1.summ, CC2.summ, res.CCP = CCP, k = 4,
                        element = "class", alg.names = paste0("A", 1:8))

z1 <- consensus_compare(x, cl.mat = y2, cons.mat = y1)
z2 <- consensus_compare(x, cl.mat = y2, cons.mat = y1,
                        alg.names =  paste0("A", 1:8))

test_that("combining results has expected lengths", {
  expect_length(y1, length(CCP) + dim(CC1)[3] + dim(CC2)[3])
  expect_equal(ncol(y2), length(CCP) + dim(CC1)[3] + dim(CC2)[3])
})

test_that("names can be overwritten", {
  expect_identical(names(y3), paste0("A", 1:8))
  expect_identical(colnames(y4), paste0("A", 1:8))
})

test_that("comparing results works", {
  expect_error(z1, NA)
  expect_error(z2, NA)
})

test_that("weighing works", {
  expect_error(consensus_weigh(z1), NA)
})

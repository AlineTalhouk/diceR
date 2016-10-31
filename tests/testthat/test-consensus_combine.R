
context("Consensus combine, evaluate, trim, and weigh")

set.seed(911)
x <- matrix(rnorm(1000), nrow = 10)
CC1 <- ConClust(x, nc = 2:4, reps = 5, method = "apEucl", progress = FALSE)
CC2 <- ConClust(x, nc = 2:4, reps = 5, method = "gmmBIC", progress = FALSE)
set.seed(1)
ref.cl <- sample(1:4, 10, replace = TRUE)

test_that("combining results has expected lengths", {
  y1 <- consensus_combine(CC1, CC2, element = "matrix")
  y2 <- consensus_combine(CC1, CC2, element = "class")
  expect_length(y1, prod(dim(CC1)[3:4]) + prod(dim(CC2)[3:4]))
  expect_equal(ncol(y2), prod(dim(CC1)[3:4]) + prod(dim(CC2)[3:4]))
})

test_that("names can be overwritten", {
  an <- paste0(rep(c("A", "B"), 3), rep(2:4, each = 2))
  y3 <- consensus_combine(CC1, CC2, element = "matrix", alg.names = an)
  y4 <- consensus_combine(CC1, CC2, element = "class", alg.names = an)
  expect_identical(names(y3), an)
  expect_identical(colnames(y4), an)
})

test_that("comparing results works", {
  z1 <- consensus_evaluate(x, k = 4, CC1, CC2, plot = FALSE)
  z2 <- consensus_evaluate(x, k = 4, CC1, CC2, ref.cl = ref.cl, plot = FALSE)
  expect_error(z1, NA)
  expect_length(z2, 3)
})

test_that("trimming (potentially) removes algorithms", {
  CC3.combined <- abind::abind(list(CC1, CC2), along = 3)
  CC3.trimmed <- consensus_trim(x, k = 4, CC1, CC2, ref.cl = ref.cl,
                                quantile = 0.8)$E_trimmed
  expect_lte(dim(CC3.trimmed)[3], dim(CC3.combined)[3])
})
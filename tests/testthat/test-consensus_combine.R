
context("Consensus combine, evaluate, trim, and weigh")

set.seed(911)
x <- matrix(rnorm(300), nrow = 100)
CC1 <- consensus_cluster(x, nk = 2:4, reps = 5, algorithms = "ap",
                         progress = FALSE)
CC2 <- consensus_cluster(x, nk = 2:4, reps = 5, algorithms = "gmm",
                         progress = FALSE)
CC3.combined <- abind::abind(list(CC1, CC2), along = 3)
ref.cl <- sample(1:4, 100, replace = TRUE)

test_that("combining results has expected lengths", {
  y1 <- consensus_combine(CC1, CC2, element = "matrix")
  y2 <- consensus_combine(CC1, CC2, element = "class")
  expect_length(unlist(y1, recursive = FALSE),
                prod(dim(CC1)[3:4]) + prod(dim(CC2)[3:4]))
  expect_equal(ncol(data.frame(y2)), prod(dim(CC1)[3:4]) + prod(dim(CC2)[3:4]))
})

test_that("evaluation works with reference class and can plot", {
  cons.cl <- matrix(sample(1:4, 400, replace = TRUE), ncol = 4,
                    dimnames = list(NULL, LETTERS[1:4]))
  expect_error(consensus_evaluate(x, CC1, CC2, cons.cl = cons.cl, plot = TRUE),
               NA)
  expect_length(consensus_evaluate(x, CC1, CC2, ref.cl = ref.cl, plot = FALSE),
                5)
})

test_that("there are different ways to choose k", {
  expect_error(consensus_evaluate(x, CC1, CC2, k.method = "all"), NA)
  expect_error(consensus_evaluate(x, CC1, CC2, k.method = 3), NA)
  expect_error(consensus_evaluate(x, CC1, CC2, k.method = 2:3))
})

test_that("compactness measure works with singleton clusters", {
  ref.cl <- c(sample(1:3, 99, replace = TRUE), 4)
  expect_error(compactness(x, ref.cl), NA)
})

test_that("trimming (potentially) removes algorithms", {
  CC3.trimmed <- consensus_evaluate(x, CC1, CC2, ref.cl = ref.cl, n = 1,
                                    trim = TRUE)$trim$data.new
  expect_lte(dim(CC3.trimmed)[3], dim(CC3.combined)[3])
})

test_that("reweighing (potentially) replicates each slice of algorithm", {
  CC3.trimmed <- consensus_evaluate(x, CC1, CC2, ref.cl = ref.cl,
                                    trim = TRUE, reweigh = TRUE)
  expect_error(CC3.trimmed, NA)
})

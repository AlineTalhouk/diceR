
context("Consensus combine, evaluate, trim, and weigh")

set.seed(911)
x <- matrix(rnorm(1000), nrow = 100)
CC1 <- consensus_cluster(x, nk = 2:4, reps = 5, algorithms = "ap", progress = FALSE)
CC2 <- consensus_cluster(x, nk = 2:4, reps = 5, algorithms = "gmm", progress = FALSE)
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
                4)
})

test_that("compactness measure works with singleton clusters", {
  ref.cl <- c(sample(1:3, 99, replace = TRUE), 4)
  expect_error(compactness(x, ref.cl), NA)
})

test_that("trimming (potentially) removes algorithms", {
  CC3.combined <- abind::abind(list(CC1, CC2), along = 3)
  CC3.trimmed <- consensus_trim(x, CC1, CC2, ref.cl = ref.cl,
                                quantile = 0.8)$data.new
  expect_lte(dim(CC3.trimmed)[3], dim(CC3.combined)[3])
})

test_that("reweighing (potentially) replicates each slice of algorithm", {
  CC3.combined <- abind::abind(list(CC1, CC2), along = 3)
  CC3.trimmed <- consensus_trim(x, CC1, CC2, ref.cl = ref.cl,
                                reweigh = TRUE)
  expect_error(CC3.trimmed, NA)
})

test_that("trimming doesn't have to show evaluation", {
  CC3.trimmed <- consensus_trim(x, CC1, CC2, ref.cl = ref.cl,
                                quantile = 0.8, show.eval = FALSE)
  expect_null(CC3.trimmed$eval)
})
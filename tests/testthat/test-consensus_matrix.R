context("Consensus matrix")

set.seed(1)
x <- replicate(100, rbinom(100, 4, 0.2))
w <- rexp(100)
w <- w / sum(w)
y1 <- consensus_matrix(x)
y2 <- consensus_matrix(x, w)

test_that("consensus matrix is square", {
  expect_equal(nrow(y1), 100)
  expect_equal(ncol(y1), 100)
})

test_that("weighting works", {
  expect_equal(nrow(y2), 100)
  expect_equal(ncol(y2), 100)
  expect_false(isTRUE(all.equal(y1, y2)))
})

test_that("missing entries are handled", {
  x[, 50] <- NA
  y3 <- consensus_matrix(x)
  expect_equal(nrow(y3), 100)
  expect_equal(ncol(y3), 100)
})

test_that("single vector computation works", {
  set.seed(1)
  sv <- rbinom(100, 4, 0.2)
  expect_error(consensus_matrix(sv), NA)
})

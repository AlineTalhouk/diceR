
context("Weigh each pair of clusters using their shared members (Jaccard coefficient)")

set.seed(1)
E <-
  matrix(rep(sample(1:4, 1000, replace = TRUE)), nrow = 100, byrow = FALSE)
mat <- matrix(c(1, 4, 3, 2, 2, 2, 1, 3, 3, 1, 4, 4, 4, 3, 2, 1), nrow = 4)
mat.wcl <- matrix(c(0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0),
                       nrow = 4)

test_that("Check weigh_clusters works with LCE", {
  expect_equal(sum(!weigh_clusters(E) == t(weigh_clusters(E))), 0)
  expect_equal(diag(weigh_clusters(E)), rep(0,nrow(weigh_clusters(E))))
  expect_equal(weigh_clusters(E)[2, 1], 0.9)
  expect_equal(weigh_clusters(E)[4, 1], 0.91)
})

test_that("It works with another matrix", {
  mat.wcl <- matrix(c(0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0),
                         nrow = 4)
  expect_equal(sum(!weigh_clusters(mat) == mat.wcl), 0)
})

test_that("Error with wrong inputs", {
  expect_error(weigh_clusters(data.frame(
    val1 = c(1, 2, 3, 4), val2 = c(5, 6, 7, 8))))
  expect_error(weigh_clusters(1))
  expect_error(weigh_clusters(c(1, 2, 3)))
  expect_error(weigh_clusters(matrix(c("a", "b", "c", "d", "e", "f"),
                                     nrow = 3)))
})

test_that("Error with negative labels", {
  mat_negLabel <- mat
  mat_negLabel[2, 1] <- -2
  expect_error(weigh_clusters(mat_negLabel))
})
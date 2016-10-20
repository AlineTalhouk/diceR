
context("Distance criterion clustering")

set.seed(1)
E <- matrix(rep(sample(1:4, 1000, replace = TRUE)), nrow = 100)
dc <- 0.8
K <- 4

test_that("Check dist_cluster with E_LCE", {
  tree_complete <- linkage(Y = stod(asrs(E, dc)), method = "complete")
  tree_complete_clustered <- dist_cluster(Z = tree_complete, maxclust = K)
  tree_average <- linkage(Y = stod(asrs(E, dc)), method = "average")
  tree_average_clustered <- dist_cluster(Z = tree_average, maxclust = K)
  tree_single <- linkage(Y = stod(asrs(E, dc)), method = "single")
  tree_single_clustered <- dist_cluster(Z = tree_single, maxclust = K)
  expect_equal(nrow(tree_complete), nrow(E) - 1)
  expect_equal(nrow(tree_average), nrow(E) - 1)
  expect_equal(nrow(tree_single), nrow(E) - 1)
  expect_equal(nrow(E), nrow(tree_average_clustered))
  expect_equal(nrow(E), nrow(tree_complete_clustered))
  expect_equal(nrow(E), nrow(tree_single_clustered))
  expect_equal(sum(sort(unique(tree_single_clustered)) != 1:K), 0)
  expect_equal(sum(sort(unique(tree_complete_clustered)) != 1:K), 0)
  expect_equal(sum(sort(unique(tree_average_clustered)) != 1:K), 0)
})

test_that("Check dist_cluster throws error with wrong inputs", {
  tree <- linkage(stod(asrs(E, dc)), "complete")
  expect_error(dist_cluster(Z = tree, maxclust = -3))
  expect_error(dist_cluster(matrix(letters[1:16], ncol = 4), 4))
  expect_error(dist_cluster(c(1, 2, 3, 4), 2))
})

test_that("Check dist_cluster with a very small tree",{
  tree_short <- linkage(stod(cts(E, dc)), "complete")[1:3, ]
  expect_equal(nrow(dist_cluster(Z = tree_short, maxclust = K)), K)
  expect_equal(sum(sort(unique(
    dist_cluster(Z = tree_short, maxclust = K))) != 1:K), 0)
  expect_equal(sum(unique(dist_cluster(Z = tree_short, maxclust = 1)) != 1), 0)
})
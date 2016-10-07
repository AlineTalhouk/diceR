
context("dist_cluster")

test_that("Check dist_cluster with E_LCE", {
  data("E_LCE")
  dc<-0.8
  K<-4
  tree<-linkage(stod(asrs(E_LCE,dc)),"complete")
  tree_clustered<-dist_cluster(Z=tree,maxclust = K)
  expect_equal(nrow(tree),nrow(E_LCE)-1)
  expect_equal(nrow(E_LCE),nrow(tree_clustered))
  expect_equal(sum(sort(unique(tree_clustered))!=1:K),0)
})

test_that("Check dist_cluster throws error with wrong inputs",{
  data("E_LCE")
  dc<-0.8
  tree<-linkage(stod(asrs(E_LCE,dc)),"complete")
  expect_error(dist_cluster(Z=tree,maxclust = -3))
  expect_error(dist_cluster(Z=tree,maxclust="three"))
  expect_error(dist_cluster(matrix(letters[1:16],ncol=4),4))
  expect_error(dist_cluster(c(1,2,3,4),2))
})
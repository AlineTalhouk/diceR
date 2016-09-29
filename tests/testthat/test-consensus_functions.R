context("Consensus functions")

test_that("majority voting works", {
  dt = array(c(2,3,2,2,2,3,1,1,2,3,2,2,2,3,1,1),c(2,4,2))
  majority_voting(dt)
  expect_equal( majority_voting(dt, is.relabelled = TRUE), c("2","3"))
  })

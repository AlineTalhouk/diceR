context("Principal component Normal procedure")

test_that("pcn simulation and selection works", {
  set.seed(9)
  A <- matrix(rnorm(300), nrow = 20)
  pc.dat <- pcn_simulate(A, n.sim = 50)
  cl <- sample(1:4, 20, replace = TRUE)
  pc.select <- pcn_select(pc.dat, cl, "rep")
  expect_length(pc.select, 2)

  pc.select <- pcn_select(pc.dat, cl, "range")
  expect_length(pc.select, 3)
})

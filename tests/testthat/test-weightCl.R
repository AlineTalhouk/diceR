
context("weight for each pair of clusters using their shared members (Jaccard coefficient)")

data(E_LCE)
mat <- matrix(c(1, 4, 3, 2, 2, 2, 1, 3, 3, 1, 4, 4, 4, 3, 2, 1), nrow =
                4)
mat.weightCl <- matrix(c(0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0), nrow =
                         4)

test_that("Check weightCl works with LCE", {

  expect_equal(sum(!weightCl(E_LCE) == t(weightCl(E_LCE))), 0)
  expect_equal(diag(weightCl(E_LCE)), c(0, 0, 0, 0))
  expect_equal(weightCl(E_LCE)[2, 1], 1)
  expect_equal(weightCl(E_LCE)[4, 1], 0.91)
})

test_that("Check wightCl with another matrix",{
  mat.weightCl <- matrix(c(0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0), nrow =
                           4)
  expect_equal(sum(!weightCl(mat) == mat.weightCl), 0)
})

test_that("Check weigthCl with wrong inputs",{
  expect_error(weightCl(data.frame(
    val1 = c(1, 2, 3, 4), val2 = c(5, 6, 7, 8)
  )))
  expect_error(weightCl(1))
  expect_error(weightCl(c(1, 2, 3)))
  expect_error(weightCl(matrix(c(
    "a", "b", "c", "d", "e", "f"
  ), nrow = 3)))
})

test_that("Check weightCl with negative labels",{
  mat_negLabel <- mat
  mat_negLabel[2, 1] <- -2
  expect_error(weightCl(mat_negLabel))
})

context("Get row column numbers of a number in a matrix")

test_that("Check getRowColNumbers",{

  mat<-matrix(c(60,17,58,62,81,11,32,7,28,85,80,15,19,50,45,40,88,31,84,30,99,94,61,55,27),ncol=5)
  expect_true(getRowColNumbers(mat,60)$rows==1)
  expect_true(getRowColNumbers(mat,60)$cols==1)
  expect_true(getRowColNumbers(mat,85)$rows==5)
  expect_true(getRowColNumbers(mat,85)$cols==2)
  expect_true(getRowColNumbers(mat,84)$rows==4)
  expect_true(getRowColNumbers(mat,84)$cols==4)
  expect_true(getRowColNumbers(mat,94)$rows==2)
  expect_true(getRowColNumbers(mat,94)$cols==5)

  mat_wd<-matrix(c(60,17,58,62,81,11,32,7,28,85,80,15,19,50,45,40,88,32,84,30,32,94,61,85,27),ncol=5)
  expect_equal(getRowColNumbers(mat_wd,32)$rows,c(2,3,1))
  expect_equal(getRowColNumbers(mat_wd,32)$cols,c(2,4,5))
  expect_equal(getRowColNumbers(mat_wd,85)$rows,c(5,4))
  expect_equal(getRowColNumbers(mat_wd,85)$cols,c(2,5))

  expect_error(getRowColNumbers(mat,101))
  expect_error(getRowColNumbers(mat,c(32,80)))
  expect_error(getRowColNumbers(mat,"101"))
  expect_error(getRowColNumbers(data.frame(val1=c(1,2,3),val2=c(4,5,6)),1))
  expect_error(getRowColNumbers(matrix(letters[1:25],ncol=5),"a"))
  expect_error(getRowColNumbers(matrix(letters[1:25],ncol=5),1))
})

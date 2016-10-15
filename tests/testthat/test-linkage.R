context("linkage")

set.seed(1)
E_LCE <-
  matrix(rep(sample(1:4, 1000, replace = TRUE)), nrow = 100, byrow = FALSE)

test_that("Check linkage with complete linkage",{
  Y<-linkage(stod(asrs(E_LCE,0.8)),"complete")
  expect_equal(nrow(Y),nrow(E_LCE)-1)
  expect_true(abs(sum(Y[,3])-36.9299)<0.0001)
})

test_that("Check linkage with average linkage",{
  Y<-linkage(stod(asrs(E_LCE,0.8)),"average")
  expect_equal(nrow(Y),nrow(E_LCE)-1)
  expect_true(abs(sum(Y[,3])-35.8230)<0.0001)
})

test_that("Check linkage with single linkage",{
  Y<-linkage(stod(asrs(E_LCE,0.8)),"single")
  expect_equal(nrow(Y),nrow(E_LCE)-1)
  expect_true(abs(sum(Y[,3])-33.1691)<0.0001)
})

test_that("Check linkage with wrong inputs",{
  expect_error(linkage(c("a","b","c"),"complete"))
  expect_error(linkage(stod(asrs(E_LCE,0.8)),"bccrc"))
  expect_error(linkage(asrs(E_LCE,0.8),"average"))
})



# test_that("Check UIJ",{
#   test_uij1<-UIJ(4,4,8)
#   expect_equal(test_uij1$U,c(1,2,3,5,6,7,8))
#   expect_equal(test_uij1$I,c(3,9,14,19,20,21,22))
#   expect_equal(test_uij1$J,c(3,9,14,19,20,21,22))
# 
#   test_uij2<-UIJ(1,4,8)
#   expect_equal(test_uij2$U,c(2,3,5,6,7,8))
#   expect_equal(test_uij2$I,c(1,2,4,5,6,7))
#   expect_equal(test_uij2$J,c(9,14,19,20,21,22))
#   
#   test_uij3<-UIJ(5,3,8)
#   expect_equal(test_uij3$U,c(1,2,3,4,4,5,6,7,8))
#   expect_equal(test_uij3$I,c(4,10,15,19,21,22,23,24,25))
#   expect_equal(test_uij3$J,c(2,8,13,17,14,15,16,17,18))
# })
# 
# test_that("Check UIJ throws error with wrong inputs",{
#   expect_error(UIJ(i=4,j=4,m="eight"))
#   expect_error(UIJ(i=4,j="four",m=8))
#   expect_error(UIJ(i="four",j=4,m=8))
#   expect_error(UIJ(i=c(1,2,3),6,8))
#   expect_error(i=4,j=c(8,7,6),m=9)
#   expect_error(i=4,j=5,m=c(11,23))
# })
context("stod, distance vector")

test_that("Check stod",{
  data(E_LCE)
  S<-cts(E_LCE,dc=0.8)
  s<-stod(S)
  expect_true(sum(!s[1:24]==0)==0)
  expect_true(abs(s[25]-0.5979)<=0.0001)
  expect_true(abs(sum(s)-2962)<1)
  expect_error(stod(1))
  expect_error(stod(c(1,2)))
  expect_error(stod("vancouver"))
  expect_error(stod(matrix(c("a","b","c","d","e","f"),nrow=3)))
})

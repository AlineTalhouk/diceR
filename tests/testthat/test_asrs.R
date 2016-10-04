test_that("Check asrs",{
  data("E_LCE")
  ASRS<-asrs(E_LCE,0.8)
  expect_equal(sum(!diag(ASRS)==1),0)
  expect_true(abs(sum(ASRS)-1456)<1)
  expect_error(asrs(E_LCE))
  expect_error(asrs(E_LCE,1.2))
  expect_error(asrs(E_LCE,-9))
  expect_error(asrs(c(1,1,2,3),0.8))

  })

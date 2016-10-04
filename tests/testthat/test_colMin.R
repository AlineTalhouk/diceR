test_that("Check colMin",{
  expect_error(colMin(data.frame(names=c("a","b","c"),vals=c(17,18,19))))
  expect_error(colMin(data.frame(states=c("CALI","TEXS","NEWY"),cities=c("san fran","houston","new york"))))
  expect_error(colMin(data.frame(val1=c(1,2,3),val2=c(4,5,6))))
  expect_equal(colMin(matrix(1:16,ncol=4)),c(1,5,9,13))
  expect_equal(colMin(matrix(c(60,17,58,62,81,11,32,7,28,85,80,15,19,50,45,40,88,31,84,30,99,94,61,55,27),
                             ncol=5)),
               c(17,7,15,30,27))
})

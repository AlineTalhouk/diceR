
context("Relabel clusters")

mat <- matrix(c(1, 2, 3, 4, 4, 2, 1, 3, 3, 1, 4, 2, 2, 3, 4, 1),
              byrow = TRUE,
              nrow = 4)

test_that("It works for relabelling and total number of classes", {
  mat_relabelled <- relabel_clusters(mat)
  expect_true(sum(!mat_relabelled$newE == matrix(
    c(1, 4, 3, 2, 6, 6, 5, 7, 9, 8, 10, 10, 14, 13, 12, 11), nrow = 4
  )) == 0)
  expect_equal(mat_relabelled$no_allcl, 14)
})

test_that("Error with negative class", {
  mat_negClass <- mat
  mat_negClass[2, 3] <- -1
  expect_error(relabel_clusters(mat_negClass))
})

test_that("Error with wrong inputs", {
  expect_error(relabel_clusters(3))
  expect_error(relabel_clusters(c(1, 2)))
  expect_error(relabel_clusters(data.frame(
    letters = c("a", "b", "c"), vals = c(1, 2, 3)
  )))
  expect_error(relabel_clusters(matrix(c("a", "b", "c"))))
})

context("External validity indices")

set.seed(1)
x <- sample(1:4, 100, replace = TRUE)
y1 <- sample(1:4, 100, replace = TRUE)
y2 <- sample(1:3, 100, replace = TRUE)

test_that("normalized mutual information works", {
  expect_error(ev_nmi(x, y1), NA)
})

test_that("error if different number of unique labels", {
  expect_error(ev_confmat(x, y1), NA)
  expect_error(ev_confmat(x, y2))
})

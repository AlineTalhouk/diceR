context("prepare_data")

data(hgsc)
hgsc <- hgsc[1:40, 1:30]

test_that("tsne and largevis options result in 2D representation", {
  dat_tsne <- prepare_data(hgsc, type = "tsne")
  dat_lv <- prepare_data(hgsc, type = "largevis")
  expect_equal(ncol(dat_tsne), 2)
  expect_equal(ncol(dat_lv), 2)
})

test_that("conventional and robust options reduce dimension but center/scale", {
  dat_con <- prepare_data(hgsc, type = "conventional")
  dat_rob <- prepare_data(hgsc, type = "robust")
  expect_lte(ncol(dat_con), ncol(hgsc))
  expect_lte(ncol(dat_rob$data), ncol(hgsc))
})

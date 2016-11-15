
context("Diverse Cluster Ensemble")

data(hgsc)
dat <- t(hgsc[, -1])

test_that("dice works with one algorithm, one consensus funs", {
  dice.obj <- dice(dat, nk = 4, algorithms = "hcAEucl", consensusFUNS = "kmodes")
  expect_length(dice.obj, 2)
  expect_equal(dim(dice.obj$clusters), c(nrow(dat), 1))
})

test_that("dice works with multiple algorithms, consensus funs, trimming, and reference class", {
  ref.cl <- data.frame(initCol = rownames(dat)) %>%
    tidyr::separate(initCol,
                    into = c("patientID", "Class"),
                    sep = "_") %>% 
    use_series(Class) %>% 
    factor() %>% 
    as.integer()
  
  dice.obj <- dice(dat, nk = 4,
                   algorithms = c("hcAEucl", "hcDianaEucl", "pamEucl", "pamSpear"),
                   consensusFUNS = c("kmodes", "majority", "LCE"),
                   trim = TRUE, ref.cl = ref.cl)
  expect_length(dice.obj, 2)
  expect_is(dice.obj$clusters, "matrix")
})
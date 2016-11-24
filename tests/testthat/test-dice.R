
context("Diverse Cluster Ensemble")

data(hgsc)
dat <- t(hgsc[, -1])

test_that("dice works with one algorithm, one consensus funs", {
  dice.obj <- dice(dat, nk = 4, algorithms = "hc", cons.funs = "kmodes")
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
  
  dice.obj <- dice(dat, nk = 4, reps = 5,
                   algorithms = c("hc", "diana", "pam"),
                   cons.funs = c("kmodes", "majority"),
                   trim = TRUE, ref.cl = ref.cl)
  expect_length(dice.obj, 2)
  expect_is(dice.obj$clusters, "matrix")
})
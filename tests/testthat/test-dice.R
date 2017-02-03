
context("Diverse Cluster Ensemble")

data(hgsc)
dat <- t(hgsc[, -1])
ref.cl <- data.frame(initCol = rownames(dat)) %>%
  tidyr::separate(initCol,
                  into = c("patientID", "Class"),
                  sep = "_") %>% 
  use_series(Class) %>% 
  factor() %>% 
  as.integer()

test_that("dice works with one algorithm, one consensus funs", {
  dice.obj <- dice(dat, nk = 4, algorithms = "hc", cons.funs = "kmodes")
  expect_length(dice.obj, 2)
  expect_equal(dim(dice.obj$clusters), c(nrow(dat), 1))
})

test_that("dice works with multiple algorithms, consensus funs, trimming, and reference class", {
  dice.obj <- dice(dat, nk = 4, reps = 5,
                   algorithms = c("hc", "diana", "pam"),
                   cons.funs = c("kmodes", "majority"),
                   trim = TRUE, ref.cl = ref.cl)
  expect_length(dice.obj, 2)
  expect_is(dice.obj$clusters, "matrix")
})

test_that("single algorithm and single consensus return same results", {
  dice.obj <- dice(dat, nk = 4, reps = 50, algorithms = "km",
                   cons.funs = "CSPA", ref.cl = ref.cl)
  ind.obj <- unname(dice.obj$indices$external)
  expect_identical(as.character(ind.obj[1, -1]), as.character(ind.obj[2, -1]))
})
context("Diverse Cluster Ensemble")

library(dplyr)
data(hgsc)
ref.cl <- strsplit(rownames(hgsc), "_") %>% 
  purrr::map_chr(2) %>% 
  factor() %>% 
  as.integer()

test_that("dice works with one algorithm, one consensus funs", {
  dice.obj <- dice(hgsc, nk = 4, algorithms = "hc", cons.funs = "kmodes")
  expect_length(dice.obj, 5)
  expect_equal(dim(dice.obj$clusters), c(nrow(hgsc), 1))
})

test_that("dice works with multiple algorithms, consensus funs, trimming, and reference class", {
  dice.obj <- dice(hgsc, nk = 4, reps = 5,
                   algorithms = c("hc", "diana", "pam"),
                   cons.funs = c("kmodes", "majority"),
                   trim = TRUE, n = 2, ref.cl = ref.cl)
  expect_length(dice.obj, 5)
  expect_is(dice.obj$clusters, "matrix")
})

test_that("single algorithm and single consensus return same results", {
  dice.obj <- dice(hgsc, nk = 4, reps = 50, algorithms = "km",
                   cons.funs = "CSPA", ref.cl = ref.cl)
  ind.obj <- unname(dice.obj$indices$external$`4`)
  expect_identical(as.character(ind.obj[1, -1]), as.character(ind.obj[2, -1]))
})

test_that("indices slot returns NULL if evaluate specified as FALSE", {
  dice.obj <- dice(hgsc, nk = 4, reps = 3, algorithms = "hc",
                   cons.funs = "kmodes", ref.cl = ref.cl, evaluate = FALSE)
  expect_null(dice.obj$indices)
})

test_that("relabelling uses 1st col if more than 1 cons.funs and no ref.cl", {
  dice.obj <- dice(hgsc, nk = 4, reps = 3, algorithms = "hc",
                   cons.funs = c("kmodes", "majority"),  evaluate = FALSE)
  expect_error(dice.obj, NA)
})

test_that("cluster size prepended when multiple k requested", {
  dice.obj <- dice(hgsc, nk = 3:4, reps = 3, algorithms = "hc",
                   cons.funs = "kmodes", k.method = "all", evaluate = FALSE)
  expect_true(all(grepl("k=", colnames(dice.obj$clusters))))
})
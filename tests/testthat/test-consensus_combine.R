skip_if_not_installed("apcluster")

set.seed(911)
x <- matrix(rnorm(300), nrow = 100)
CC1 <- consensus_cluster(x, nk = 2:4, reps = 5, algorithms = "ap",
                         progress = FALSE)
CC2 <- consensus_cluster(x, nk = 2:4, reps = 5, algorithms = "gmm",
                         progress = FALSE)
CC3 <- consensus_cluster(x, nk = 2:4, reps = 5, algorithms = "hc",
                         progress = FALSE)
ref.cl <- sample(1:4, 100, replace = TRUE)

consensus_summary_reference <- function(E) {
  hc_fn <- getFromNamespace("hc", "diceR")
  con.mats <- E %>%
    purrr::array_tree(c(4, 3)) %>%
    purrr::modify_depth(2, consensus_matrix) %>%
    purrr::map(magrittr::set_names, dimnames(E)[[3]]) %>%
    magrittr::set_names(dimnames(E)[[4]])
  con.cls <- con.mats %>%
    purrr::imap(~ purrr::map(.x, function(z) hc_fn(stats::dist(z), k = .y)))
  dplyr::lst(con.mats, con.cls) %>% purrr::transpose()
}

consensus_combine_reference <- function(..., element = c("matrix", "class")) {
  cs <- abind::abind(list(...), along = 3) %>%
    consensus_summary_reference()
  switch(
    match.arg(element),
    matrix = purrr::map(cs, "con.mats"),
    class = purrr::map(cs, "con.cls") %>%
      purrr::map(~ do.call(cbind, .))
  )
}

test_that("combining results has expected lengths", {
  y1 <- consensus_combine(CC1, CC2, element = "matrix")
  y2 <- consensus_combine(CC1, CC2, element = "class")
  expect_length(unlist(y1, recursive = FALSE),
                prod(dim(CC1)[3:4]) + prod(dim(CC2)[3:4]))
  expect_equal(ncol(data.frame(y2)), prod(dim(CC1)[3:4]) + prod(dim(CC2)[3:4]))
})

test_that("evaluation works with reference class and can plot", {
  cons.cl <- matrix(sample(1:4, 400, replace = TRUE), ncol = 4,
                    dimnames = list(NULL, LETTERS[1:4]))
  expect_length(consensus_evaluate(x, CC1, CC2, cons.cl = cons.cl,
                                   ref.cl = ref.cl, plot = TRUE),
                5)
})

test_that("there are different ways to choose k", {
  expect_error(consensus_evaluate(x, CC1, CC2, k.method = "all"), NA)
  expect_error(consensus_evaluate(x, CC1, CC2, k.method = 3), NA)
  expect_error(consensus_evaluate(x, CC1, CC2, k.method = 2:3))
})

test_that("compactness measure works with singleton clusters", {
  ref.cl <- c(sample(1:3, 99, replace = TRUE), 4)
  expect_error(compactness(x, ref.cl), NA)
})

test_that("trimming (potentially) removes algorithms", {
  CC.trimmed <- consensus_evaluate(x, CC1, CC2, ref.cl = ref.cl, n = 1,
                                   trim = TRUE)$trim.obj$E.new
  expect_lte(dim(CC.trimmed[[1]])[3],
             dim(abind::abind(list(CC1, CC2), along = 3))[3])
})

test_that("reweighing (potentially) replicates each slice of algorithm", {
  CC.trimmed1 <- consensus_evaluate(x, CC1, CC2, ref.cl = ref.cl,
                                    trim = TRUE, reweigh = TRUE,
                                    k.method = "all")
  CC.trimmed2 <- consensus_evaluate(x, CC1, CC2, CC3, ref.cl = ref.cl,
                                    trim = TRUE, reweigh = TRUE, n = 2)
  expect_error(CC.trimmed1, NA)
  expect_error(CC.trimmed2, NA)
})

test_that("consensus_combine lazy class path matches reference behavior", {
  args <- list(CC1, CC2)
  ref_mat <- do.call(consensus_combine_reference, c(args, list(element = "matrix")))
  new_mat <- do.call(consensus_combine, c(args, list(element = "matrix")))
  expect_equal(new_mat, ref_mat, tolerance = 1e-12)

  ref_cls <- do.call(consensus_combine_reference, c(args, list(element = "class")))
  new_cls <- do.call(consensus_combine, c(args, list(element = "class")))
  expect_equal(new_cls, ref_cls, tolerance = 1e-12)
})

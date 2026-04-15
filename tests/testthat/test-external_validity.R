set.seed(1)
x <- sample(1:4, 100, replace = TRUE)
y1 <- sample(1:4, 100, replace = TRUE)
y2 <- sample(1:3, 100, replace = TRUE)

concordance_counts_reference <- function(cl1, cl2) {
  n <- length(cl1)
  yy <- yn <- ny <- nn <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      if (cl1[[i]] == cl1[[j]]) {
        if (cl2[[i]] == cl2[[j]]) {
          yy <- yy + 1
        } else {
          yn <- yn + 1
        }
      } else {
        if (cl2[[i]] == cl2[[j]]) {
          ny <- ny + 1
        } else {
          nn <- nn + 1
        }
      }
    }
  }
  list(yy = yy, yn = yn, ny = ny, nn = nn)
}

test_that("normalized mutual information works", {
  expect_error(ev_nmi(x, y1), NA)
})

test_that("error if different number of unique labels", {
  expect_error(ev_confmat(x, y1), NA)
  expect_error(ev_confmat(x, y2))
})

test_that("concordance counts closed-form matches reference loops", {
  set.seed(9)
  cl1 <- sample(1:7, 300, replace = TRUE)
  cl2 <- sample(1:5, 300, replace = TRUE)
  cc <- getFromNamespace("concordance_counts", "diceR")
  expect_identical(cc(cl1, cl2), concordance_counts_reference(cl1, cl2))
})

test_that("concordance counts rejects missing labels", {
  cc <- getFromNamespace("concordance_counts", "diceR")
  expect_error(
    cc(c(1L, 1L, NA_integer_, 2L), c(1L, 2L, 1L, 2L)),
    "must not contain missing values"
  )
})

test_that("pair-count external indices match reference concordance", {
  set.seed(15)
  cl1 <- sample(1:6, 250, replace = TRUE)
  cl2 <- sample(1:4, 250, replace = TRUE)
  cc <- concordance_counts_reference(cl1, cl2)
  Nt <- Reduce(`+`, cc)

  hubert_ref <- (Nt * cc[["yy"]] - (cc[["yy"]] + cc[["yn"]]) *
                   (cc[["yy"]] + cc[["ny"]])) /
    sqrt((cc[["yy"]] + cc[["yn"]]) * (cc[["yy"]] + cc[["ny"]]) *
           (cc[["nn"]] + cc[["yn"]]) * (cc[["nn"]] + cc[["ny"]]))
  jaccard_ref <- cc[["yy"]] / (cc[["yy"]] + cc[["yn"]] + cc[["ny"]])
  mcnemar_ref <- (cc[["yn"]] - cc[["ny"]]) / sqrt(cc[["yn"]] + cc[["ny"]])
  rand_ref <- (cc[["yy"]] + cc[["nn"]]) / Nt

  ev_hubert_fn <- getFromNamespace("ev_hubert", "diceR")
  ev_jaccard_fn <- getFromNamespace("ev_jaccard", "diceR")
  ev_mcnemar_fn <- getFromNamespace("ev_mcnemar", "diceR")
  ev_rand_fn <- getFromNamespace("ev_rand", "diceR")

  expect_equal(ev_hubert_fn(cl1, cl2), hubert_ref, tolerance = 1e-12)
  expect_equal(ev_jaccard_fn(cl1, cl2), jaccard_ref, tolerance = 1e-12)
  expect_equal(ev_mcnemar_fn(cl1, cl2), mcnemar_ref, tolerance = 1e-12)
  expect_equal(ev_rand_fn(cl1, cl2), rand_ref, tolerance = 1e-12)
})

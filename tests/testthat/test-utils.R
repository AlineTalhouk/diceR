cl <- seq_len(4)

min_fnorm_reference <- function(A, B = diag(nrow(A))) {
  n <- nrow(A)
  D <- matrix(NA, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      D[j, i] <- sum((B[j, ] - A[i, ]) ^ 2)
    }
  }
  vec <- clue::solve_LSAP(D)
  list(pmat = A[vec, ], perm = vec, ord = order(vec))
}

test_that("relabelling outputs a integer", {
  set.seed(2)
  pred <- sample(cl, 100, replace = TRUE)
  true <- sample(cl, 100, replace = TRUE)
  expect_type(relabel_class(pred, true), "integer")
})

test_that("flatten uses first clustering as reference if not relabelled", {
  E <- matrix(sample(cl, 500, replace = TRUE), ncol = 5)
  E4 <- array(sample(cl, 2000, replace = TRUE), dim = c(100, 5, 2, 2))
  expect_error(flatten_E(E, is.relabelled = FALSE), NA)
  expect_error(flatten_E(E4, is.relabelled = FALSE), NA)
})

test_that("Check is_pos_int works", {
  expect_true(is_pos_int(3))
  expect_false(is_pos_int(-3))
  expect_true(is_pos_int(1e6))
  expect_false(is_pos_int(3.6))
  expect_false(is_pos_int(3.21e1))
  expect_false(is_pos_int(-3.7))
  expect_false(is_pos_int(-3.21e1))
  expect_false(is_pos_int(0))
})

test_that("Error if input is not a number or more than one element", {
  expect_error(is_pos_int(c(1, 2, 3)))
  expect_error(is_pos_int("a"))
  expect_error(is_pos_int(matrix(c(1, 2, 3), nrow = 1)))
})

test_that("min_fnorm matrix algebra matches reference behavior", {
  set.seed(101)
  A <- matrix(rnorm(49), nrow = 7)
  B <- matrix(rnorm(49), nrow = 7)
  ref_default <- min_fnorm_reference(A)
  new_default <- min_fnorm(A)
  expect_equal(new_default$pmat, ref_default$pmat, tolerance = 1e-12)
  expect_equal(new_default$perm, ref_default$perm, tolerance = 1e-12)
  expect_identical(new_default$ord, ref_default$ord)

  ref_custom <- min_fnorm_reference(A, B)
  new_custom <- min_fnorm(A, B)
  expect_equal(new_custom$pmat, ref_custom$pmat, tolerance = 1e-12)
  expect_equal(new_custom$perm, ref_custom$perm, tolerance = 1e-12)
  expect_identical(new_custom$ord, ref_custom$ord)
})

test_that("min_fnorm handles numerically delicate self-matches", {
  set.seed(1)
  for (iter in 1:3) {
    A <- matrix(1e4 + rnorm(6), nrow = 2)
  }
  a2 <- rowSums(A ^ 2)
  raw_D <- matrix(a2, nrow = nrow(A), ncol = nrow(A)) +
    matrix(a2, nrow = nrow(A), ncol = nrow(A), byrow = TRUE) -
    2 * (A %*% t(A))

  expect_lte(abs(min(raw_D)), 1e-6)

  ref <- min_fnorm_reference(A, A)
  expect_no_error(new <- min_fnorm(A, A))
  expect_equal(new$pmat, ref$pmat, tolerance = 1e-12)
  expect_equal(new$perm, ref$perm, tolerance = 1e-12)
  expect_identical(new$ord, ref$ord)
})

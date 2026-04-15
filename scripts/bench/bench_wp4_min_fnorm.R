#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  load_dicer <- function() {
    if ("package:diceR" %in% search()) return(invisible(TRUE))
    if (requireNamespace("pkgload", quietly = TRUE)) {
      loaded <- try(pkgload::load_all(".", export_all = FALSE, quiet = TRUE),
                    silent = TRUE)
      if (!inherits(loaded, "try-error")) return(invisible(TRUE))
    }
    if (requireNamespace("diceR", quietly = TRUE)) {
      library(diceR)
      return(invisible(TRUE))
    }
    stop("Could not load diceR for benchmarking.")
  }
  load_dicer()
})

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

elapsed_times <- function(fun, n_iter) {
  out <- numeric(n_iter)
  for (i in seq_len(n_iter)) {
    gc(verbose = FALSE)
    out[i] <- unname(system.time(fun())[["elapsed"]])
  }
  out
}

set.seed(20260210)
A <- matrix(rnorm(250 * 250), nrow = 250)
B <- matrix(rnorm(250 * 250), nrow = 250)

ref_fun <- function() min_fnorm_reference(A)
new_fun <- function() min_fnorm(A)

warmup_ref_default <- min_fnorm_reference(A)
warmup_new_default <- min_fnorm(A)
stopifnot(isTRUE(all.equal(warmup_new_default, warmup_ref_default,
                           tolerance = 1e-12)))
stopifnot(isTRUE(all.equal(min_fnorm(A, B), min_fnorm_reference(A, B),
                           tolerance = 1e-12)))

n_iter <- 7L
time_ref <- elapsed_times(ref_fun, n_iter = n_iter)
time_new <- elapsed_times(new_fun, n_iter = n_iter)

median_ref <- stats::median(time_ref)
median_new <- stats::median(time_new)
speedup <- median_ref / median_new

cat("WP4 Benchmark: min_fnorm matrix algebra rewrite (default B path)\n")
cat("Runs:", n_iter, "\n")
cat(sprintf("Baseline median (s): %.6f\n", median_ref))
cat(sprintf("New median (s): %.6f\n", median_new))
cat(sprintf("Speedup (baseline/new): %.3fx\n", speedup))

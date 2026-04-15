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

elapsed_times <- function(fun, n_iter) {
  out <- numeric(n_iter)
  for (i in seq_len(n_iter)) {
    gc(verbose = FALSE)
    out[i] <- unname(system.time(fun())[["elapsed"]])
  }
  out
}

set.seed(20260210)
cl1 <- sample(1:12, 5000, replace = TRUE)
cl2 <- sample(1:9, 5000, replace = TRUE)
cc_new <- getFromNamespace("concordance_counts", "diceR")

ref_fun <- function() concordance_counts_reference(cl1, cl2)
new_fun <- function() cc_new(cl1, cl2)

warmup_ref <- ref_fun()
warmup_new <- new_fun()
stopifnot(isTRUE(all.equal(warmup_new, warmup_ref, tolerance = 1e-12)))

n_iter <- 7L
time_ref <- elapsed_times(ref_fun, n_iter = n_iter)
time_new <- elapsed_times(new_fun, n_iter = n_iter)

median_ref <- stats::median(time_ref)
median_new <- stats::median(time_new)
speedup <- median_ref / median_new

cat("WP5 Benchmark: concordance_counts closed-form counts\n")
cat("Runs:", n_iter, "\n")
cat(sprintf("Baseline median (s): %.6f\n", median_ref))
cat(sprintf("New median (s): %.6f\n", median_new))
cat(sprintf("Speedup (baseline/new): %.3fx\n", speedup))

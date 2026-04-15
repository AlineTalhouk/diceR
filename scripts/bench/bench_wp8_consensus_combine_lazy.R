#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (requireNamespace("magrittr", quietly = TRUE)) {
    library(magrittr)
  }
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

elapsed_times <- function(fun, n_iter) {
  out <- numeric(n_iter)
  for (i in seq_len(n_iter)) {
    gc(verbose = FALSE)
    out[i] <- unname(system.time(fun())[["elapsed"]])
  }
  out
}

set.seed(20260210)
x <- matrix(rnorm(220 * 60), nrow = 220, ncol = 60)
cc1 <- consensus_cluster(x, nk = 2:4, reps = 8, algorithms = "hc", progress = FALSE)
cc2 <- consensus_cluster(x, nk = 2:4, reps = 8, algorithms = "pam", progress = FALSE)

ref_fun <- function() consensus_combine_reference(cc1, cc2, element = "matrix")
new_fun <- function() consensus_combine(cc1, cc2, element = "matrix")

warmup_ref <- ref_fun()
warmup_new <- new_fun()
stopifnot(isTRUE(all.equal(warmup_new, warmup_ref, tolerance = 1e-12)))

n_iter <- 7L
time_ref <- elapsed_times(ref_fun, n_iter = n_iter)
time_new <- elapsed_times(new_fun, n_iter = n_iter)

median_ref <- stats::median(time_ref)
median_new <- stats::median(time_new)
speedup <- median_ref / median_new

cat("WP8 Benchmark: consensus_combine lazy class computation\n")
cat("Runs:", n_iter, "\n")
cat(sprintf("Baseline median (s): %.6f\n", median_ref))
cat(sprintf("New median (s): %.6f\n", median_new))
cat(sprintf("Speedup (baseline/new): %.3fx\n", speedup))

#' Parallel consensus function
#'
#' The consensus functions are passed from `dice`.
#'
#' @inheritParams dice
#'
#' @return A vector of final consensus class assignments using one of the
#'   consensus functions.
#' @export
pl_confun <- function(nk, cons.funs = c("kmodes", "majority", "CSPA", "LCE"),
                      sim.mat = c("cts", "srs", "asrs")) {
  Ecomp <- readRDS(fs::dir_ls("merged", regexp = "merged/Ecomp.rds"))
  E <- readRDS(fs::dir_ls("merged", regexp = "merged/E.rds"))

  fs::dir_create("consensus")
  purrr::walk(cons.funs, ~ {
    purrr::map(sim.mat, function(sm) {
      consensus <- switch(
        .,
        kmodes = k_modes(Ecomp),
        majority = majority_voting(Ecomp),
        CSPA = CSPA(E, nk),
        LCE = LCE(Ecomp, k = nk, sim.mat = sm)
      )
      fn <- ifelse(. == "LCE", paste0(., sm), .)
      saveRDS(consensus, paste0("consensus/cons_", fn, ".rds"))
    })
  })
}

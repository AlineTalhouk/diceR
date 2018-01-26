#' Parallel cluster merging
#'
#' Produce the objects E, Eknn, and Ecomp from `dice` and save to file.
#'
#' @inheritParams consensus_cluster
#'
#' @return A raw merged cluster array, a knn-imputed cluster array, and a
#'   completely imputed (using majority voting) cluster array all saved to file
#' @export
pl_merge <- function(data, nk = 4) {
  E <- pl_merge_impl("raw")
  Eknn <- pl_merge_impl("imputed")
  Ecomp <- impute_missing(Eknn, data, nk = nk)
  fs::dir_create("merged")
  dplyr::lst(E, Eknn, Ecomp) %>%
    purrr::set_names(paste0("merged/", names(.), ".rds")) %>%
    purrr::iwalk(saveRDS)
}

#' Parallel merge helper
#' @noRd
pl_merge_impl <- function(path) {
  out <- fs::dir_ls(path) %>%
    purrr::map(readRDS) %>%
    split(purrr::map_chr(purrr::map(., dimnames), 3)) %>%
    purrr::map(abind::abind, along = 2) %>%
    abind::abind(along = 3)
  dimnames(out)[[2]] <- paste0("R", seq_along(dimnames(out)[[2]]))
  out
}

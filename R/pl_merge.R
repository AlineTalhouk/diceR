#' Parallel cluster and consensus matrix merging
#'
#' `pl_merge_clust` constructs the objects E, Eknn, and Ecomp from `dice` and
#' saves to file. `pl_merge_conmat` constructs consensus matrices for each
#' algorithm and saves to file.
#'
#' @inheritParams consensus_cluster
#'
#' @return `pl_merge_clust` returns a raw merged cluster array, a knn-imputed
#'   cluster array, and a completely imputed (using majority voting) cluster
#'   array all saved to file. `pl_merge_conmat` returns a list of consensus
#'   matrices.
#' @name pl_merge
#' @export
pl_merge_clust <- function(data, nk = 4) {
  E <- pl_merge_clust_impl("raw")
  Eknn <- pl_merge_clust_impl("imputed")
  Ecomp <- impute_missing(Eknn, data, nk = nk)
  fs::dir_create("merged")
  dplyr::lst(E, Eknn, Ecomp) %>%
    purrr::set_names(paste0("merged/", names(.), ".rds")) %>%
    purrr::iwalk(saveRDS)
}

#' Parallel merge cluster helper
#' @noRd
pl_merge_clust_impl <- function(path) {
  out <- fs::dir_ls(path) %>%
    purrr::map(readRDS) %>%
    split(purrr::map_chr(purrr::map(., dimnames), 3)) %>%
    purrr::map(abind::abind, along = 2) %>%
    abind::abind(along = 3)
  dimnames(out)[[2]] <- paste0("R", seq_along(dimnames(out)[[2]]))
  out
}

#' @rdname pl_merge
#' @export
pl_merge_conmat <- function() {
  fs::dir_create("merged")
  conmat <- fs::dir_ls("conmat") %>%
    purrr::map(readRDS) %>% {
      lapply(rapply(., enquote, how = "unlist"), eval)
    } %>%
    split(stringr::str_split_fixed(names(.), "\\.", n = 4)[, 4]) %>%
    lapply(pl_merge_conmat_impl)
  meta_conmat <- pl_merge_conmat_impl(conmat)
  dplyr::lst(conmat, meta_conmat) %>%
    purrr::set_names(paste0("merged/", names(.), ".rds")) %>%
    purrr::iwalk(saveRDS)
}

#' Parallel merge consensus matrices helper
#' @noRd
pl_merge_conmat_impl <- function(x) {
  Reduce("+", purrr::map(x, ~ . / length(x)))
}

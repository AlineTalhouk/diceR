#' Parallel K-Nearest Neighbours imputation
#'
#' Impute each individual clustering unit from `pl_cluster` and write results to
#' a new directory.
#'
#' @inheritParams impute_knn
#' @return one file for each "raw" output object, this time prefixed with
#'   "imputed".
#' @export
pl_impute <- function(data, seed = 123) {
  fs::dir_create("imputed")
  fs::dir_ls("raw") %>%
    purrr::map(readRDS) %>%
    purrr::map(apply, 2:4, impute_knn, data = data, seed = seed) %>%
    purrr::set_names(gsub("raw", "imputed", names(.))) %>%
    purrr::iwalk(saveRDS)
}

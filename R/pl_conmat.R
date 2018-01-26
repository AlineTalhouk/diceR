#' Parallel consensus matrix
#'
#' Create a consensus matrix for each imputed cluster vector from the parallel
#' pipeline.
#'
#' @return one file for each "imputed" output object, this time prefixed with
#'   "conmat".
#' @export
pl_conmat <- function() {
  fs::dir_create("conmat")
  fs::dir_ls("imputed") %>%
    purrr::map(readRDS) %>%
    purrr::map(consensus_combine, element = "matrix") %>%
    purrr::set_names(gsub("imputed", "conmat", names(.))) %>%
    purrr::iwalk(saveRDS)
}

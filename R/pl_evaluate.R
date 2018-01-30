#' Parallel evaluation of clustering results
#'
#' Internal evaluation indices, cluster assignments, and summaries are returned.
#'
#' The `CC` is the most commonly selected cluster. The `Pr` is the proportion at
#' which the cluster `CC` appears across different algorithms.
#'
#' @inheritParams dice
#' @return Three objects: The internal evaluation indices matrix, a data frame
#'   of cluster assignments for each of the algorithms and consensus functions,
#'   and a data frame of the top `n` selected algorithms plus the `CC` and `Pr`
#'   summaries.
#' @export
pl_evaluate <- function(data, nk, n = 5) {
  # Read in the consensus matrices
  cons.mat <- fs::dir_ls("conmat") %>%
    purrr::map(readRDS)
  cl.mat <- readRDS("merged/conmat.rds") %>%
    purrr::map(~ hc(stats::dist(.), k = nk)) %>%
    purrr::map_df(relabel_class, ref.cl = .[[1]])
  final <- fs::dir_ls("consensus") %>%
    purrr::map_df(readRDS) %>%
    dplyr::mutate_all(as.integer) %>%
    purrr::set_names(gsub(".*cons_(.*).rds", "\\1", names(.))) %>%
    cbind(cl.mat)

  # relabel the elements of the data frame
  finalR <- final
  finalR[] <- apply(final, 2, relabel_class, ref.cl = final[, "majority"])

  # Cluster evaluate at this point
  ii <- ivi_table(finalR, data)

  # Separate algorithms into those from clusterCrit (main), and (others)
  cr <- consensus_rank(ii, n = n)
  top <- cr$top.list
  ii <- ii[match(top, ii$Algorithms), ]
  finalR <- finalR[, top]
  saveRDS(ii, "evaluate/ii.rds")
  saveRDS(finalR, "evaluate/all_clusts.rds")

  # Add CC and Pr summaries to Final output
  Final <- finalR %>%
    dplyr::select(top[seq_len(n)]) %>%
    cbind(apply(., 1, function(x) {
      CC <- names(which.max(table(x)))
      Pr <- table(x)[CC] / length(x)
      list(CC = as.integer(CC), Pr = Pr)
    }) %>%
      dplyr::bind_rows())
  saveRDS(Final, "evaluate/FinalR.rds")
}

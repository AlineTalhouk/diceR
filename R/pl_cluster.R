#' Parallel consensus_cluster
#'
#' @inheritParams consensus_cluster
#' @param ... other arguments passed to `consensus_cluster`
#'
#' @return one file for each combination of `algorithms`, `distance`, and
#'   `reps`. For NMF, each `nmf.method` gets a separate file. Each object is a
#'   4D array with one column, so that they can be run in parallel and combined
#'   at the end.
#' @export
#'
#' @examples
#' \dontrun{
#' data(hgsc)
#' pl_cluster(data, nk = 4, reps = 2, algorithms = c("nmf", "km", "block"),
#' distance = "manhattan", nmf.method = "lee", file.name = "test")
#' }
pl_cluster <- function(data, nk = 4, reps = 1000, algorithms = NULL,
                       nmf.method = c("brunet", "lee"), distance = "euclidean",
                       file.name = NULL, ...) {
  args <- dplyr::lst(data, nk, reps = 1, ...)
  r <- seq_len(reps)
  algs <- dplyr::lst(NALG, DALG, OALG) %>%
    purrr::map(~ algorithms[algorithms %in% .])
  fs::dir_create("raw")
  pl_cluster_nmf(args, r, algs$NALG, nmf.method, file.name)
  pl_cluster_dist(args, r, algs$DALG, distance, file.name)
  pl_cluster_other(args, r, algs$OALG, file.name)
}

#' Parallel clustering for NMF algorithms
#' @noRd
pl_cluster_nmf <- function(args, reps, algorithms, nmf.method, file.name) {
  if (length(algorithms) > 0) {
    purrr::map(nmf.method, function(m) {
      purrr::map(reps, function(r) {
        purrr::invoke(
          consensus_cluster,
          args,
          algorithms = algorithms,
          nmf.method = m,
          seed.data = r,
          file.name = paste0("raw/",
                             paste("raw", algorithms, m, r, file.name, sep = "_"))
        )
      })
    })
  }
}

#' Parallel clustering for distance-based algorithms
#' @noRd
pl_cluster_dist <- function(args, reps, algorithms, distance, file.name) {
  if (length(algorithms) > 0) {
    purrr::map(algorithms, function(a) {
      purrr::map(distance, function(d) {
        purrr::map(reps, function(r) {
          purrr::invoke(
            consensus_cluster,
            args,
            algorithms = a,
            seed.data = r,
            file.name = paste0("raw/",
                               paste("raw", a, d, r, file.name, sep = "_"))
          )
        })
      })
    })
  }
}

#' Parallel clustering for other algorithms
#' @noRd
pl_cluster_other <- function(args, reps, algorithms, file.name) {
  if (length(algorithms) > 0) {
    purrr::map(algorithms, function(a) {
      purrr::map(reps, function(r) {
        purrr::invoke(
          consensus_cluster,
          args,
          algorithms = a,
          seed.data = r,
          file.name = paste0("raw/",
                             paste("raw", a, r, file.name, sep = "_"))
        )
      })
    })
  }
}

#' Diverse Clustering Ensemble
#'
#' Runs consensus clustering across subsamples, algorithms, and number of
#' clusters (k).
#'
#' There are three ways to handle the input data before clustering via argument
#' \code{prep.data}. The default is to use the raw data as-is ("none"). Or, we
#' can enact \code{\link{prepare_data}} on the full dataset ("full"), or the
#' bootstrap sampled datasets ("sampled").
#'
#' @param cons.funs consensus functions to use. Current options are "kmodes"
#'   (k-modes), "majority" (majority voting), "CSPA" (Cluster-based Similarity
#'   Partitioning Algorithm), "LCE" (linkage clustering ensemble)
#' @param evaluate logical; if \code{TRUE} (default), validity indices are
#'   returned. Internal validity indices are always computed. If \code{ref.cl}
#'   is not \code{NULL}, then external validity indices will also be computed.
#' @inheritParams consensus_cluster
#' @inheritParams consensus_evaluate
#' @inheritParams LCE
#' @inheritParams impute_knn
#' @return A list with the following elements
#' \item{E}{raw clustering ensemble object}
#' \item{Eknn}{clustering ensemble object with knn imputation used on \code{E}}
#' \item{Ecomp}{flattened ensemble object with remaining missing entries imputed
#' by majority voting}
#' \item{clusters}{final clustering assignment from the diverse clustering
#' ensemble method}
#' \item{indices}{if \code{evaluate = TRUE}, shows cluster evaluation indices;
#' otherwise \code{NULL}}
#' @author Aline Talhouk, Derek Chiu
#' @export
#' @examples
#' library(dplyr)
#' data(hgsc)
#' dat <- hgsc[1:100, 1:50]
#' ref.cl <- strsplit(rownames(dat), "_") %>%
#'   purrr::map_chr(2) %>%
#'   factor() %>%
#'   as.integer()
#' dice.obj <- dice(dat, nk = 4, reps = 5, algorithms = "hc", cons.funs =
#' "kmodes", ref.cl = ref.cl, progress = FALSE)
#' str(dice.obj, max.level = 2)
dice <- function(data, nk, reps = 10, algorithms = NULL, k.method = NULL,
                 nmf.method = c("brunet", "lee"), distance = "euclidean",
                 cons.funs = c("kmodes", "majority", "CSPA", "LCE"),
                 sim.mat = c("cts", "srs", "asrs"),
                 prep.data = c("none", "full", "sampled"), min.var = 1,
                 seed = 1, trim = FALSE, reweigh = FALSE, n = 5,
                 evaluate = TRUE, plot = FALSE, ref.cl = NULL,
                 progress = TRUE) {

  # Check that inputs are correct
  assertthat::assert_that(length(dim(data)) == 2)
  prep.data <- match.arg(prep.data)

  # Generate Diverse Cluster Ensemble
  E <- consensus_cluster(data = data, nk = nk, reps = reps,
                         algorithms = algorithms, nmf.method = nmf.method,
                         distance = distance, prep.data = prep.data,
                         min.var = min.var, progress = progress)
  # KNN imputation
  Eknn <- apply(E, 2:4, impute_knn, data = data, seed = seed)

  # Select k and new (trimmed and reweighed) data
  if (progress)
    cli::cat_line("Selecting k and imputing non-clustered cases")
  eval.obj <- consensus_evaluate(data = data, Eknn, ref.cl = ref.cl,
                                 k.method = k.method, trim = trim,
                                 reweigh = reweigh, n = n)
  Eknn <- eval.obj$trim$data.new
  k <- eval.obj$k

  # Impute remaining missing cases
  Ecomp <- purrr::map2(Eknn, k, impute_missing, data = data)

  # Consensus functions
  if (progress)
    cli::cat_line("Computing consensus functions")
  Final <- purrr::map2(Ecomp, k, ~ {
    vapply(cons.funs, function(x) {
      switch(x,
             kmodes = k_modes(.x),
             majority = majority_voting(.x),
             CSPA = CSPA(E, .y),
             LCE = LCE(drop(.x), k = .y, sim.mat = sim.mat)
      )
    }, double(nrow(.x))) %>%
      apply(2, as.integer)
    })

  #  If more than one k, need to prepend "k=" labels
  if (length(Ecomp) > 1) {
    Final <- purrr::map2(Final, k,
                         ~ magrittr::set_colnames(.x, paste0(colnames(.),
                                                             " k=", .y)))
  }

  # Relabel Final Clustering using reference (or first column if no reference)
  if (is.null(ref.cl)) {
    FinalR <- purrr::map(Final, ~ apply(.x, 2, relabel_class, ref.cl = .x[, 1]))
  } else {
    FinalR <- purrr::map(Final, ~ apply(.x, 2, relabel_class, ref.cl = ref.cl))
  }
  clusters <- FinalR %>%
    purrr::invoke(cbind, .) %>%
    magrittr::set_rownames(rownames(data))

  # Return evaluation output including consensus function results
  if (evaluate) {
    if (progress)
      cli::cat_line("Evaluating output with consensus function results")
    eval.obj2 <- consensus_evaluate(data, E, cons.cl = clusters,
                                    ref.cl = ref.cl, plot = plot)
    indices <- c(k = list(eval.obj[["k"]]), eval.obj2[2:4],
                 trim = list(eval.obj[["trim"]]))
  } else {
    indices <- NULL
  }

  # Add the reference class as the first column if provided
  if (!is.null(ref.cl)) {
    clusters <- cbind(Reference = ref.cl, clusters)
  }

  # Remove list structure
  Eknn <- abind::abind(Eknn, along = 3)
  Ecomp <- abind::abind(Ecomp, along = 3)
  if (progress)
    cli::cat_line("Diverse Cluster Ensemble Completed")
  dplyr::lst(E, Eknn, Ecomp, clusters, indices)
}

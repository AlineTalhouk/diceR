#' Diverse Clustering Ensemble
#'
#' Runs consensus clustering across subsamples, algorithms, and number of
#' clusters (k).
#'
#' There are three ways to handle the input data before clustering via argument
#' `prep.data`. The default is to use the raw data as-is ("none"). Or, we can
#' enact [prepare_data()] on the full dataset ("full"), or the bootstrap sampled
#' datasets ("sampled").
#'
#' @param cons.funs consensus functions to use. Current options are "kmodes"
#'   (k-modes), "majority" (majority voting), "CSPA" (Cluster-based Similarity
#'   Partitioning Algorithm), "LCE" (linkage clustering ensemble), "LCA" (latent
#'   class analysis)
#' @param evaluate logical; if `TRUE` (default), validity indices are returned.
#'   Internal validity indices are always computed. If `ref.cl` is not `NULL`,
#'   then external validity indices will also be computed.
#' @param plot logical; if `TRUE`, `graph_all` is called and a summary
#'   evaluation heatmap of ranked algorithms vs. internal validity indices is
#'   plotted as well.
#' @inheritParams consensus_cluster
#' @inheritParams consensus_evaluate
#' @inheritParams LCE
#' @inheritParams impute_knn
#' @return A list with the following elements
#' \item{E}{raw clustering ensemble object}
#' \item{Eknn}{clustering ensemble object with knn imputation used on `E`}
#' \item{Ecomp}{flattened ensemble object with remaining missing entries imputed
#' by majority voting}
#' \item{clusters}{final clustering assignment from the diverse clustering
#' ensemble method}
#' \item{indices}{if `evaluate = TRUE`, shows cluster evaluation indices;
#' otherwise `NULL`}
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
dice <- function(data, nk, p.item = 0.8, reps = 10, algorithms = NULL, k.method = NULL,
                 nmf.method = c("brunet", "lee"), hc.method = "average",
                 distance = "euclidean",
                 cons.funs = c("kmodes", "majority", "CSPA", "LCE", "LCA"),
                 sim.mat = c("cts", "srs", "asrs"),
                 prep.data = c("none", "full", "sampled"), min.var = 1,
                 seed = 1, seed.data = 1, trim = FALSE, reweigh = FALSE, n = 5,
                 evaluate = TRUE, plot = FALSE, ref.cl = NULL,
                 progress = TRUE) {

  # Check that inputs are correct
  assertthat::assert_that(length(dim(data)) == 2)
  prep.data <- match.arg(prep.data)

  # Generate Diverse Cluster Ensemble
  E <- consensus_cluster(data = data, nk = nk, p.item = p.item, reps = reps,
                         algorithms = algorithms, nmf.method = nmf.method,
                         hc.method = hc.method, distance = distance,
                         prep.data = prep.data, min.var = min.var,
                         seed.data = seed.data, progress = progress)
  # KNN imputation
  Eknn <- apply(E, 2:4, impute_knn, data = data, seed = seed)

  # Select k and new (trimmed and reweighed) data
  if (progress)
    cat("Selecting k and imputing non-clustered cases\n")

  eval.obj <- consensus_evaluate(data = data, Eknn, ref.cl = ref.cl,
                                 k.method = k.method, trim = trim,
                                 reweigh = reweigh, n = n)
  trim.obj <- eval.obj$trim.obj
  Eknn <- trim.obj$E.new
  k <- eval.obj$k

  # Impute remaining missing cases
  Ecomp <- purrr::map2(Eknn, k, impute_missing, data = data)

  # Consensus functions
  if (progress)
    cat("Computing consensus functions\n")

  Final <- purrr::map2(Ecomp, k, ~ {
    vapply(cons.funs, function(x) {
      switch(x,
             kmodes = k_modes(.x),
             majority = majority_voting(.x),
             CSPA = CSPA(E, .y),
             LCE = LCE(drop(.x), k = .y, sim.mat = sim.mat),
             LCA = LCA(.x)
      ) %>% as.integer()
    }, integer(nrow(.x)))
  })

  #  If more than one k, need to prepend "k=" labels
  if (length(Ecomp) > 1)
    Final <- purrr::map2(Final, k,
                         ~ magrittr::set_colnames(.x, paste_k(colnames(.), .y)))

  # Relabel Final Clustering using reference (or first column if no reference)
  clusters <- Final %>%
    purrr::map(~ apply(., 2, relabel_class, ref.cl = ref.cl %||% .[, 1])) %>%
    do.call(cbind, .) %>%
    magrittr::set_rownames(rownames(data))

  # Return evaluation output including consensus function results
  if (evaluate) {
    if (progress)
      cat("Evaluating output with consensus function results\n")

    eval.obj2 <- consensus_evaluate(data = data, E, cons.cl = clusters,
                                    ref.cl = ref.cl, plot = plot)
    indices <- c(k = k, eval.obj2[2:4], trim = list(trim.obj))
  } else {
    indices <- NULL
  }

  # Add the reference class as the first column if provided
  if (!is.null(ref.cl))
    clusters <- cbind(Reference = ref.cl, clusters)

  # Algorithm vs internal index heatmap
  if (plot)
    algii_heatmap(data, k, E, clusters, ref.cl)

  # Combine DICE with different E and indices
  if (progress)
    cat("Diverse Cluster Ensemble Completed\n")

  dplyr::lst(E, Eknn, Ecomp, clusters, indices) %>%
    purrr::map_at(2:3, abind::abind, along = 3) # Unlist Eknn/Ecomp
}

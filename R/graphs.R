#' Graphical Displays
#'
#' Graph cumulative distribution function (CDF) graphs, relative change in area
#' under CDF curves, heatmaps, and cluster assignment tracking plots.
#'
#' \code{graph_cdf} plots the CDF for consensus matrices from different
#' algorithms. \code{graph_delta_area} calculates the relative change in area
#' under CDF curve between algorithms. \code{graph_heatmap} generates consensus
#' matrix heatmaps for each algorithm in \code{x}. \code{graph_tracking} tracks
#' how cluster assignments change between algorithms. \code{graph_all} is a
#' wrapper that runs all graphing functions.
#'
#' @param x an object from \code{\link{consensus_cluster}}
#' @param mat same as \code{x}, or a list of consensus matrices computed from
#'   \code{x} for faster results
#' @param cl same as \code{x}, or a matrix of consensus classes computed from
#'   \code{x} for faster results
#' @return Various plots from \code{graph_*{}} functions. All plots are
#'   generated using \code{ggplot}, except for \code{graph_heatmap}, which uses
#'   \code{\link[gplots]{heatmap.2}}. Colours used in \code{graph_heatmap} and
#'   \code{graph_tracking} utilize \code{\link[RColorBrewer]{RColorBrewer}}
#'   palettes.
#' @name graphs
#' @author Derek Chiu
#' @export
#' @examples
#' # Consensus clustering for 3 algorithms
#' library(ggplot2)
#' set.seed(911)
#' x <- matrix(rnorm(100), ncol = 10)
#' CC1 <- consensus_cluster(x, nk = 2:4, reps = 5,
#' algorithms = c("hc", "ap", "gmm"), progress = FALSE)
#'
#' # Plot CDF
#' p <- graph_cdf(CC1)
#'
#' # Change y label and add colours
#' p + labs(y = "Probability") + stat_ecdf(aes(colour = k)) +
#' scale_color_brewer(palette = "Set2")
#'
#' # Delta Area
#' p <- graph_delta_area(CC1)
#'
#' # Heatmaps with column side colours corresponding to clusters
#' CC2 <- consensus_cluster(x, nk = 3, reps = 3, algorithms = "hc", progress =
#' FALSE)
#' graph_heatmap(CC2)
#'
#' # Track how cluster assignments change between algorithms
#' p <- graph_tracking(CC1)
graph_cdf <- function(mat) {
  dat <- get_cdf(mat)
  p <- ggplot(dat, aes_(x = ~CDF, colour = ~k)) +
    stat_ecdf() +
    facet_wrap(~Method) +
    labs(x = "Consensus Index",
         y = "CDF",
         title = "Consensus Cumulative Distribution Functions")
  print(p)
  return(p)
}

#' @rdname graphs
#' @export
graph_delta_area <- function(mat) {
  dat <- get_cdf(mat) %>%
    dplyr::group_by_("Method", "k") %>%
    dplyr::summarize_(.dots = stats::setNames(
      list(~ flux::auc(seq(0, 1, length.out = table(k)[1]), CDF)), "AUC")) %>%
    dplyr::mutate_(.dots = stats::setNames(
      list(~c(AUC[1], diff(AUC) / AUC[-length(AUC)])), "da"))
  if (length(unique(dat$k)) > 1) {
    p <- ggplot(dat, aes_(x = ~k, y = ~da)) +
      geom_line(group = 1) +
      geom_point() +
      facet_wrap(~Method) +
      labs(y = "Relative change in Area under CDF curve",
           title = "Delta Area")
    print(p)
    return(p)
  }
}

#' Calculate CDF for each clustering algorithm at each k
#' @noRd
get_cdf <- function(mat) {
  if (inherits(mat, "array")) {
    mat <- consensus_combine(mat, element = "matrix")
  }
  mat %>%
    purrr::modify_depth(2, ~ .x[lower.tri(.x, diag = TRUE)]) %>%
    purrr::imap(~ purrr::set_names(.x, paste(.y, names(.x), sep = "."))) %>%
    dplyr::bind_cols() %>%
    tidyr::gather_("Group", "CDF", names(.)) %>%
    tidyr::separate_("Group", c("k", "Method"), sep = "\\.")
}

#' @param main heatmap title. If \code{NULL} (default), the titles will be taken
#'   from names in \code{mat}
#' @param ... additional arguments to \code{\link[gplots]{heatmap.2}}
#'
#' @rdname graphs
#' @export
graph_heatmap <- function(mat, main = NULL, ...) {
  if (inherits(mat, "array")) {
    mat <- consensus_combine(mat, element = "matrix")
  }
  dat <- mat %>%
    purrr::flatten() %>%
    magrittr::set_names(list(purrr::map(mat, names)[[1]], names(mat)) %>%
                          purrr::cross() %>%
                          purrr::map_chr(paste, collapse = " k="))
  main <- paste(main %||% names(dat), "Consensus Matrix")
  assertthat::assert_that(length(main) == length(purrr::flatten(mat)))

  cs.col <- RColorBrewer::brewer.pal(8, "Set2")
  cc <- purrr::map2(dat, rep(as.numeric(names(mat)),
                             each = unique(purrr::map_int(mat, length))),
                    ~ sample(cs.col)[hc(stats::dist(.x), k = .y)])
  purrr::pwalk(list(dat, main, cc), function(dat, main, cc)
    gplots::heatmap.2(x = dat, main = main, ColSideColors = cc,
                      col = grDevices::colorRampPalette(
                        RColorBrewer::brewer.pal(n = 9, "PuBuGn"))(256),
                      labRow = "", labCol = "", trace = "none",
                      hclustfun = function(d)
                        stats::hclust(d, method = "average"),
                      dendrogram = "column", ...))
}

#' @rdname graphs
#' @export
graph_tracking <- function(cl) {
  if (inherits(cl, "array")) {
    cl <- consensus_combine(cl, element = "class")
  }
  dat <- cl %>%
    purrr::imap(~ `colnames<-`(.x, paste(.y, colnames(.x), sep = "."))) %>%
    do.call(cbind, .) %>%
    as.data.frame() %>%
    tidyr::gather_("Group", "Class", names(.)) %>%
    tidyr::separate_("Group", c("k", "Method"), sep = "\\.") %>%
    cbind(Samples = seq_len(unique(purrr::map_int(cl, nrow)))) %>%
    dplyr::mutate_at(dplyr::vars(c("Class", "Method", "Samples")), factor)
  if (length(unique(dat$k)) > 1) {
    p <- ggplot(dat, aes_(x = ~Samples, y = ~k)) +
      geom_tile(aes_(fill = ~Class)) +
      facet_wrap(~Method) +
      scale_fill_brewer(palette = "Set2") +
      ggtitle("Tracking Cluster Assignments Across k") +
      theme(axis.text.x = element_blank(),
            axis.ticks.x = element_blank())
    print(p)
    return(p)
  }
}

#' @rdname graphs
#' @export
graph_all <- function(x, ...) {
  mat <- consensus_combine(x, element = "matrix")
  cl <- consensus_combine(x, element = "class")
  graph_cdf(mat)
  graph_delta_area(mat)
  graph_heatmap(mat, ...)
  graph_tracking(cl)
}

#' Comparing ranked Algorithms vs internal indices (ii) in heatmap
#' @inheritParams dice
#' @param E object in \code{dice}
#' @param clusters object in \code{dice}
#' @noRd
algii_heatmap <- function(data, nk, E, clusters, ref.cl = NULL) {
  # Final clusters
  fc <- E %>%
    consensus_combine(element = "class") %>%
    magrittr::extract2(as.character(nk)) %>%
    cbind.data.frame(clusters) %>%
    purrr::map_df(relabel_class, ref.cl = ref.cl %||% .[, 1])

  # Internal indices
  ii <- data.frame(
    Algorithms = colnames(fc),
    fc %>% purrr::map_df(
      clusterCrit::intCriteria,
      traj = as.matrix(data),
      crit = c("Calinski_Harabasz", "Dunn", "PBM", "Tau", "Gamma", "C_index",
               "Davies_Bouldin", "McClain_Rao", "SD_Dis", "Ray_Turi", "G_plus",
               "Silhouette", "S_Dbw")),
    Compactness = fc %>% purrr::map_dbl(compactness, data = data),
    Connectivity = fc %>% purrr::map_dbl(
      ~ clValid::connectivity(Data = data, clusters = .))) %>%
    dplyr::mutate_all(dplyr::funs(structure(., names = colnames(fc))))

  # Heatmap data: reorder rows by ranked algorithms, remove indices with NaN
  hm <- ii %>%
    tibble::column_to_rownames("Algorithms") %>%
    magrittr::extract(match(rank_alg(ii), rownames(.)),
                      purrr::map_lgl(., ~ all(!is.nan(.x))))

  # Plot heatmap with annotated colours, column scaling, no further reordering
  NMF::aheatmap(
    hm,
    annCol = data.frame(Criteria = c(rep("Maximized", 5), rep("Minimized", 8))),
    annColors = list(Criteria = stats::setNames(c("darkgreen", "deeppink4"),
                                                c("Maximized", "Minimized"))),
    Colv = NA, Rowv = NA, scale = "column", col = "PiYG",
    main = "Ranked Algorithms on Internal Validity Indices"
  )
}

#' Rank Algorithms
#' @noRd
rank_alg <- function(ii) {
  # Separate internal indices into those from clusterCrit and from others
  ii.cc <- ii %>%
    magrittr::extract(!names(.) %in% c("Algorithms", "Compactness",
                                       "Connectivity") &
                        purrr::map_lgl(., ~ all(!is.nan(.x)))) # Remove NaN idx
  ii.other <- ii[c("Compactness", "Connectivity")]

  # Which algorithm is the best for each index?
  bests <- purrr::imap_int(ii.cc, clusterCrit::bestCriterion)
  max.bests <- ii.cc %>%
    magrittr::extract(purrr::map_int(., which.max) == bests) %>%
    magrittr::multiply_by(-1)
  min.bests <- ii.cc %>%
    magrittr::extract(purrr::map_int(., which.min) == bests) %>%
    cbind(ii.other)

  # Determine trimmed ensemble using rank aggregation
  rank.agg <- cbind(max.bests, min.bests) %>%
    scale(center = FALSE, scale = TRUE) %>%
    as.data.frame() %>%
    purrr::map_df(~ ii$Algorithms[order(.x, sample(length(.x)))]) %>%
    t()
  RankAggreg::RankAggreg(rank.agg, ncol(rank.agg),
                         method = "GA", verbose = FALSE)$top.list
}

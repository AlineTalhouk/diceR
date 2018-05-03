#' Graphical Displays
#'
#' Graph cumulative distribution function (CDF) graphs, relative change in area
#' under CDF curves, heatmaps, and cluster assignment tracking plots.
#'
#' `graph_cdf` plots the CDF for consensus matrices from different algorithms.
#' `graph_delta_area` calculates the relative change in area under CDF curve
#' between algorithms. `graph_heatmap` generates consensus matrix heatmaps for
#' each algorithm in `x`. `graph_tracking` tracks how cluster assignments change
#' between algorithms. `graph_all` is a wrapper that runs all graphing
#' functions.
#'
#' @param x an object from [consensus_cluster()]
#' @param mat same as `x`, or a list of consensus matrices computed from `x` for
#'   faster results
#' @param cl same as `x`, or a matrix of consensus classes computed from `x` for
#'   faster results
#' @return Various plots from \code{graph_*{}} functions. All plots are
#'   generated using `ggplot`, except for `graph_heatmap`, which uses
#'   [gplots::heatmap.2()]. Colours used in `graph_heatmap` and `graph_tracking`
#'   utilize [RColorBrewer::RColorBrewer()] palettes.
#' @name graphs
#' @author Derek Chiu
#' @export
#' @examples
#' # Consensus clustering for 3 algorithms
#' library(ggplot2)
#' set.seed(911)
#' x <- matrix(rnorm(80), ncol = 10)
#' CC1 <- consensus_cluster(x, nk = 2:4, reps = 3,
#' algorithms = c("hc", "pam", "km"), progress = FALSE)
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

#' @param main heatmap title. If `NULL` (default), the titles will be taken from
#'   names in `mat`
#' @param ... additional arguments to [gplots::heatmap.2()]
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
#' @param E object in `dice`
#' @param clusters object in `dice`
#' @noRd
algii_heatmap <- function(data, nk, E, clusters, ref.cl = NULL) {
  # Cluster list to keep
  cl.list <- E %>%
    consensus_combine(element = "class") %>%
    magrittr::extract(as.character(nk))

  # Final cluster object construction depends on value of nk
  if (length(nk) > 1) {
    fc <- purrr::map2(cl.list, nk,
                      ~ magrittr::set_colnames(.x, paste_k(colnames(.), .y))) %>%
      purrr::map2(split_clusters(clusters), cbind) %>%
      purrr::map(~ apply(., 2, relabel_class, ref.cl = ref.cl %||% .[, 1])) %>%
      do.call(cbind, .) %>%
      as.data.frame()
  } else {
    fc <- cl.list %>%
      do.call(cbind, .) %>%
      cbind.data.frame(clusters) %>%
      purrr::map_df(relabel_class, ref.cl = ref.cl %||% .[, 1])
  }

  # Internal indices
  ii <- ivi_table(fc, data)

  # Heatmap: order algorithms by ranked ii, remove indices with NaN
  hm <- ii %>%
    tibble::column_to_rownames("Algorithms") %>%
    magrittr::extract(match(consensus_rank(ii, 1)$top.list, rownames(.)),
                      purrr::map_lgl(., ~ all(!is.nan(.x))))

  # Plot heatmap with annotated colours, column scaling, no further reordering
  NMF::aheatmap(
    hm,
    annCol = data.frame(Criteria = c(rep("Maximized", 5),
                                     rep("Minimized", ncol(hm) - 5))),
    annColors = list(Criteria = stats::setNames(c("darkgreen", "deeppink4"),
                                                c("Maximized", "Minimized"))),
    Colv = NA, Rowv = NA, scale = "column", col = "PiYG",
    main = "Ranked Algorithms on Internal Validity Indices"
  )
}

#' Split clusters matrix into list based on value of k
#' @noRd
split_clusters <- function(clusters) {
  tc <- t(clusters)
  split.data.frame(x = tc,
                   f = stringr::str_split_fixed(
                     string = rownames(tc),
                     pattern = " ",
                     n = 2
                   )[, 2]) %>%
    purrr::map(t)
}

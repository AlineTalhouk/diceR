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
#' CC2 <- consensus_cluster(x, nk = 3, reps = 5, algorithms = "ap", progress =
#' FALSE)
#' graph_heatmap(CC2)
#' 
#' # Track how cluster assignments change between algorithms
#' p <- graph_tracking(CC1)
graph_cdf <- function(mat) {
  k <- CDF <- NULL
  dat <- get_cdf(mat)
  p <- ggplot(dat, aes(x = CDF, colour = k)) +
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
  k <- CDF <- AUC <- da <- Method <- NULL
  dat <- get_cdf(mat) %>% 
    group_by(Method, k) %>% 
    summarize(AUC = flux::auc(seq(0, 1, length.out = table(k)[1]), CDF)) %>% 
    mutate(da = c(AUC[1], diff(AUC) / AUC[-length(AUC)]))
  if (length(unique(dat$k)) > 1) {
    p <- ggplot(dat, aes(k, da)) +
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
  k <- Group <- CDF <- NULL
  if (inherits(mat, "array")) {
    mat <- consensus_combine(mat, element = "matrix")
  }
  dat <- mat %>% 
    purrr::at_depth(2, ~ .x[lower.tri(.x, diag = TRUE)]) %>% 
    as.data.frame() %>% 
    tidyr::gather(key = Group, value = CDF, everything()) %>%
    tidyr::separate(Group, c("k", "Method"), sep = "\\.") %>%
    mutate(k = substring(k, first = 2))
  return(dat)
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
    set_names(list(purrr::map(mat, names)[[1]], names(mat)) %>% 
                purrr::cross_n() %>% 
                purrr::map_chr(paste, collapse = " k="))
  if (is.null(main)) {
    main <- names(dat)
  } else {
    assertthat::assert_that(length(main) == length(purrr::flatten(mat)))
  }
  hm.col <- grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(n = 9, "PuBuGn"))(256)
  cs.col <- RColorBrewer::brewer.pal(8, "Set2")
  cc <- purrr::map2(dat, rep(as.numeric(names(mat)),
                             each = unique(purrr::map_int(mat, length))),
                    ~ sample(cs.col)[hc(stats::dist(.x), k = .y)])
  purrr::pwalk(list(dat, main, cc), function(dat, main, cc)
    gplots::heatmap.2(x = dat, main = paste(main, "Consensus Matrix"),
                      hclustfun = function(d)
                        stats::hclust(d, method = "average"),
                      trace = "none", dendrogram = "column", col = hm.col,
                      labRow = "", labCol = "", ColSideColors = cc, ...))
}

#' @rdname graphs
#' @export
graph_tracking <- function(cl) {
  k <- Group <- Method <- Class <- Samples <- NULL
  if (inherits(cl, "array")) {
    cl <- consensus_combine(cl, element = "class")
  }
  dat <- cl %>% 
    as.data.frame() %>% 
    tidyr::gather(key = Group, value = Class, dplyr::everything()) %>% 
    tidyr::separate(Group, c("k", "Method"), sep = "\\.") %>%
    mutate(k = substring(k, first = 2),
           Class = factor(Class), Method = factor(Method)) %>% 
    cbind(Samples = factor(seq_len(unique(purrr::map_int(cl, nrow)))),
                           levels = seq_len(unique(purrr::map_int(cl, nrow))))
  if (length(unique(dat$k)) > 1) {
    p <- ggplot(dat, aes(Samples, k)) +
      geom_tile(aes(fill = Class)) +
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
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
#' @param x an object from \code{\link{ConClust}}, or a list of multiple objects
#'   from \code{\link{ConClust}}
#' @return Various plots from \code{graph_*{}} functions. All plots are 
#'   generated using \code{ggplot}, except for \code{graph_heatmap}, which uses 
#'   \code{\link[gplots]{heatmap.2}}. Colours used in \code{graph_heatmap} and 
#'   \code{graph_tracking} utilize \code{\link[RColorBrewer]{RColorBrewer}}
#'   palettes.
#' @note The current implementations of \code{graph_delta_area} and 
#'   \code{graph_tracking} do not compare between different number of clusters. 
#'   Thus the interpretation of relative change is not natural, and the tracking
#'   plot shows low concordance.
#' @name graphs
#' @author Derek Chiu
#' @export
#' @examples
#' # Consensus clustering for 3 algorithms
#' library(ggplot2)
#' set.seed(911)
#' x <- matrix(rnorm(100), ncol = 10)
#' CC1 <- ConClust(x, nc = 2:4, reps = 5,
#' method = c("hcAEucl", "apEucl", "gmmBIC"), progress = FALSE)
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
#' CC2 <- ConClust(x, nc = 3, reps = 5, method = "apEucl", progress = FALSE)
#' graph_heatmap(CC2)
#' 
#' # Track how cluster assignments change between algorithms
#' p <- graph_tracking(CC1)
graph_cdf <- function(x) {
  k <- CDF <- NULL
  dat <- get_cdf(x)
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
graph_delta_area <- function(x) {
  k <- CDF <- AUC <- da <- Method <- NULL
  dat <- get_cdf(x) %>% 
    group_by(Method, k) %>% 
    summarize(AUC = flux::auc(seq(0, 1, length.out = table(k)[1]), CDF)) %>% 
    mutate(da = c(AUC[1], diff(AUC) / AUC[-length(AUC)]))
  p <- ggplot(dat, aes(k, da)) +
    geom_line(group = 1) +
    geom_point() +
    facet_wrap(~Method) +
    labs(y = "Relative change in Area under CDF curve",
         title = "Delta Area")
  print(p)
  return(p)
}

#' Calculate CDF for each clustering algorithm at each k
#' @noRd
get_cdf <- function(x) {
  Group <- CDF <- NULL
  assertthat::assert_that(inherits(x, "array"))
  dat <- consensus_combine(x, element = "matrix") %>% 
    lapply(function(k) lapply(k, function(d) d[lower.tri(d, diag = TRUE)])) %>% 
    as.data.frame() %>% 
    tidyr::gather(key = Group, value = CDF, dplyr::everything()) %>% 
    tidyr::separate(Group, c("k", "Method"), sep = "\\.") %>% 
    mutate(k = substring(k, first = 2))
  return(dat)
}

#' @param main heatmap title. If \code{NULL} (default), the titles will be taken
#'   from names in \code{x}
#' @param ... additional arguments to \code{\link[gplots]{heatmap.2}}
#' 
#' @rdname graphs
#' @export
graph_heatmap <- function(x, main = NULL, ...) {
  assertthat::assert_that(inherits(x, "array"))
  dat <- consensus_combine(x, element = "matrix") %>% 
    unlist(recursive = FALSE) %>% 
    set_names(names(.) %>% 
                stringr::str_split_fixed("\\.", n = 2) %>% 
                apply(1, function(x) paste0(x[2], " k=", x[1])))
  if (is.null(main)) {
    main <- names(dat)
  } else {
    assertthat::assert_that(length(main) == dim(x)[3] * dim(x)[4])
  }
  hm.col <- grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(n = 9, "PuBuGn"))(256)
  cc <- mapply(function(d, k) RColorBrewer::brewer.pal(8, "Set2")[
    stats::cutree(stats::hclust(stats::dist(d), method = "average"), k = k)],
    d = dat, k = rep(as.numeric(dimnames(x)[[4]]), each = dim(x)[3]),
    SIMPLIFY = FALSE)
  hm <- mapply(function(dat, main, cc)
    gplots::heatmap.2(x = dat, main = paste(main, "Consensus Matrix"),
                      hclustfun = function(x) stats::hclust(x, method = "average"),
                      trace = "none", dendrogram = "column", col = hm.col,
                      labRow = "", labCol = "", ColSideColors = cc, ...),
    dat = dat, main = main, cc = cc)
  return(hm)
}

#' @rdname graphs
#' @export
graph_tracking <- function(x) {
  k <- Group <- Method <- Class <- Samples <- NULL
  assertthat::assert_that(inherits(x, "array"))
  dat <- consensus_combine(x, element = "class") %>% 
    as.data.frame() %>% 
    tidyr::gather(key = Group, value = Class, dplyr::everything()) %>% 
    tidyr::separate(Group, c("k", "Method"), sep = "\\.") %>%
    mutate(k = substring(k, first = 2),
           Class = factor(Class), Method = factor(Method)) %>% 
    cbind(Samples = factor(colnames(x[, , 1, 1]),
                           levels = colnames(x[, , 1, 1])))
  p <- ggplot(dat, aes(Samples, k)) +
    geom_tile(aes(fill = Class)) +
    facet_wrap(~Method) +
    scale_fill_brewer(palette = "Set2") +
    ggtitle("Tracking Cluster Assignments Across Algorithms") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  print(p)
  return(p)
}

#' @rdname graphs
#' @export
graph_all <- function(x, ...) {
  graph_cdf(x)
  graph_delta_area(x)
  graph_heatmap(x, ...)
  graph_tracking(x)
}
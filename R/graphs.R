#' Graphical Displays
#' 
#' Graph cumulative distribution function (CDF) graphs, relative change in area 
#' under CDF curves, heatmaps, and cluster assignment tracking plots.
#' 
#' \code{graph_cdf} plots the CDF for consensus matrices from different 
#' algorithms. \code{graph_delta_area} calculates the relative change in area 
#' under CDF curve between algorithms. \code{graph_heatmap} generates consensus 
#' matrix heatmaps for each algorithm in \code{x}. \code{graph_tracking} tracks 
#' how cluster assignments change between algorithms.
#' 
#' @param x an object from \code{\link{ConClust}}
#'   
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
# Consensus clustering for 3 algorithms
#' library(ggplot2)
#' set.seed(911)
#' x <- matrix(rnorm(1000), ncol = 10)
#' CC1 <- ConClust(x, k = 4, reps = 10,
#' method = c("hcAEucl", "apEucl", "gmmBIC"))
#' 
#' # Plot CDF
#' p <- graph_cdf(CC1)
#' p
#' 
#' # Change y label and add colours
#' p + labs(y = "Probability") + stat_ecdf(aes(colour = Method)) +
#' scale_color_brewer(palette = "Set2")
#' 
#' # Delta Area
#' p <- graph_delta_area(CC1)
#' p
#' 
#' # Heatmaps with column side colours corresponding to clusters
#' CC2 <- ConClust(x, k = 3, reps = 5, method = "apEucl")
#' graph_heatmap(CC2)
#' 
#' # Track how cluster assignments change between algorithms
#' p <- graph_tracking(CC1)
#' p
graph_cdf <- function(x) {
  CDF <- NULL
  dat <- get_cdf(x)
  p <- ggplot(dat, aes(x = CDF)) +
    stat_ecdf() +
    facet_wrap(~Method) +
    labs(x = "Consensus Index",
         y = "CDF",
         title = "Consensus Cumulative Distribution Function")
  return(p)
}

#' @rdname graphs
#' @export
graph_delta_area <- function(x) {
  Method <- CDF <- AUC <- da <- NULL
  dat <- get_cdf(x)
  xvec <- seq(0, 1, length.out = table(dat$Method)[1])
  auc <- dat %>% 
    group_by(Method) %>% 
    summarize(AUC = flux::auc(xvec, CDF)) %>% 
    use_series(AUC)
  deltaK <- data.frame(Method = unique(dat$Method),
                       da = c(auc[1], diff(auc) /auc[-length(auc)]))
  p <- ggplot(deltaK, aes(Method, da)) +
    geom_line(group = 1) +
    geom_point() +
    labs(y = "Relative change in Area under CDF curve",
         title = "Delta Area")
  return(p)
}

#' Calculate CDF for each clustering algorithm
#' @noRd
get_cdf <- function(x) {
  Method <- CDF <- NULL
  dat <- consensus_summary(x, k = n_distinct(x[, 1, 1], na.rm = TRUE),
                           progress = FALSE) %>% 
    consensus_combine(element = "matrix") %>% 
    lapply(function(d) d[lower.tri(d, diag = TRUE)]) %>% 
    as.data.frame() %>% 
    tidyr::gather(key = Method, value = CDF, dplyr::everything())
  return(dat)
}

#' @param main heatmap title. If \code{NULL} (default), the titles will be taken
#'   from names in \code{x}
#' @param ... additional arguments to \code{\link[gplots]{heatmap.2}}
#' 
#' @rdname graphs
#' @export
graph_heatmap <- function(x, main = NULL, ...) {
  k <- n_distinct(x[, 1, 1], na.rm = TRUE)
  dat <- consensus_summary(x, k = k, progress = FALSE) %>% 
    consensus_combine(element = "matrix") 
  if (is.null(main)) {
    main <- names(dat)
  } else {
    assertthat::assert_that(length(main) == dim(x)[3])
  }
  hm.col <- grDevices::colorRampPalette(
    RColorBrewer::brewer.pal(n = 9, "PuBuGn"))(256)
  cc <- lapply(dat, function(d) RColorBrewer::brewer.pal(k, "Set2")[
    cutree(hclust(dist(d), method = "average"), k = k)])
  hm <- mapply(function(dat, main, cc)
    gplots::heatmap.2(x = dat, main = paste(main, "Consensus Matrix"),
                      hclustfun = function(x) hclust(x, method = "average"),
                      trace = "none", dendrogram = "column", col = hm.col,
                      labRow = "", labCol = "", ColSideColors = cc, ...),
    dat = dat, main = main, cc = cc)
  return(hm)
}

#' @rdname graphs
#' @export
graph_tracking <- function(x) {
  Method <- Class <- Samples <- NULL
  k <- n_distinct(x[, 1, 1], na.rm = TRUE)
  dat <- consensus_summary(x, k = k, progress = FALSE) %>% 
    consensus_combine(element = "class") %>% 
    as.data.frame() %>% 
    tidyr::gather(key = Method, value = Class, dplyr::everything()) %>% 
    mutate(Class = factor(Class), Method = factor(Method)) %>% 
    cbind(Samples = factor(colnames(x[, , 1]), levels = colnames(x[, , 1])))
  p <- ggplot(dat, aes(Samples, Method)) +
    geom_tile(aes(fill = Class)) +
    scale_fill_brewer(palette = "Set2") +
    ggtitle("Tracking Cluster Assignments Across Algorithms") +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())
  return(p)
}
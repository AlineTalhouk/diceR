#' Cumulative distribution function graphs
#' 
#' Graph CDF for consensus matrices from different algorithms
#' 
#' @param x an object from \code{\link{ConClust}}
#'   
#' @return a plot showing the cdfs for different algortihms
#' @export
#' @author Derek Chiu
#'   
#' @examples
#' # Consensus clustering for 3 algorithms
#' set.seed(911)
#' x <- matrix(rnorm(1000), nrow = 10)
#' CC1 <- ConClust(x, k = 4, reps = 10,
#' method = c("hcAEucl", "apEucl", "gmmBIC"), save = FALSE)
#' 
#' # Plot CDF
#' p <- graph_cdf(CC1)
#' p
#' 
#' # Change y label and add colours
#' p + labs(y = "Probability") + stat_ecdf(aes(colour = Method)) +
#' scale_color_brewer(palette = "Set2")
#' 
#' # Heatmaps with column side colours corresponding to clusters for each algorithm
#' graph_heatmap(CC1)
graph_cdf <- function(x) {
  dat <- consensus_summary(x, k = n_distinct(x[, 1, 1], na.rm = TRUE)) %>% 
    consensus_combine(element = "matrix") %>% 
    lapply(function(d) d[lower.tri(d, diag = TRUE)]) %>% 
    as.data.frame() %>% 
    tidyr::gather(key = Method, value = CDF, dplyr::everything())
  p <- ggplot(dat, aes(x = CDF)) +
    stat_ecdf() +
    facet_wrap(~Method) +
    labs(x = "Consensus Index",
         y = "CDF",
         title = "Consensus Cumulative Distribution Function")
  return(p)
}

#' Consensus matrix heatmaps
#' 
#' @param x an object from \code{\link{ConClust}}
#' @param main heatmap title
#' @param ... additional arguments to \code{\link[gplots]{heatmap.2}}
#' 
#' @return consensus matrix heatmaps for each algorithm in \code{x}
#' @author Derek Chiu
#' @export
graph_heatmap <- function(x, main = NULL, ...) {
  k <- n_distinct(x[, 1, 1], na.rm = TRUE)
  dat <- consensus_summary(x, k = k) %>% 
    consensus_combine(element = "matrix") 
  if (is.null(main)) {
    main <- names(dat)
  } else {
    assertthat::assert_that(length(main) == dim(x)[3])
  }
  hm.col <- colorRampPalette(brewer.pal(n = 9, "PuBuGn"))(256)
  cc <- lapply(dat, function(d) brewer.pal(k, "Set2")[cutree(as.dendrogram(
    hclust(dist(d), method = "average")), k = k)])
  hm <- mapply(function(dat, main, cc)
    gplots::heatmap.2(x = dat, main = paste(main, "Consensus Matrix"),
                      hclustfun = function(x) hclust(x, method = "average"),
                      trace = "none", dendrogram = "column", col = hm.col,
                      labRow = "", labCol = "", ColSideColors = cc, ...),
    dat = dat, main = main, cc = cc)
  return(hm)
}

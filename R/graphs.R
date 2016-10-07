#' Cumulative distribution function graphs
#' 
#' Graph CDF for consensus matrices from different algorithms
#' 
#' @param x an object from \code{consensus_combine} using option \code{element
#'   = "matrix"}. A list of consensus matrices from different clustering
#'   algorithms is the input.
#'   
#' @return a ggplot showing the cdfs for different algortihms
#' @export
#' @author Derek Chiu
#'   
#' @examples
#' # Consensus clustering for 3 algorithms
#' set.seed(911)
#' x <- matrix(rnorm(1000), nrow = 10)
#' CC1 <- ConClust(x, k = 4, reps = 10,
#' method = c("hcAEucl", "apEucl", "gmmBIC"), save = FALSE)
#' CC1.summ <- consensus_summary(CC1, k = 4)
#' y <- consensus_combine(CC1.summ, element = "matrix")
#' 
#' # Plot CDF
#' p <- graph_cdf(y)
#' p
#' 
#' # Change y label and add colours
#' p + labs(y = "Probability") + stat_ecdf(aes(colour = Method)) +
#' scale_color_brewer(palette = "Set2")
graph_cdf <- function(x) {
  dat <- lapply(x, function(d) d[lower.tri(d, diag = TRUE)]) %>% 
    as.data.frame() %>% 
    tidyr::gather(key = Method, value = CDF, dplyr::everything())
  p <- ggplot(dat, aes(x = CDF)) +
    stat_ecdf() +
    facet_wrap(~Method) +
    labs(y = "Prob",
         title = "Cumulative Distribution Function for each Consensus Matrix")
  return(p)
}
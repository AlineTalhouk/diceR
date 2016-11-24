#' Consensus summary
#' 
#' Given an object from \code{\link{consensus_cluster}}, returns a list of
#' consensus matrices and consensus classes for each clustering algorithm.
#' 
#' @param E output from \code{consensus_cluster}
#' @return A list with summaries for each algorithm. Each algorithm has a list 
#'   with two elements: consensus_matrix and consensus_class
#' @author Derek Chiu
#' @export
#' @examples
#' data(hgsc)
#' dat <- t(hgsc[, -1])
#' x <- consensus_cluster(dat, nk = 3:4, reps = 5, algorithms = "hc",
#' progress = FALSE)
#' cs <- consensus_summary(x)
#' str(cs)
consensus_summary <- function(E) {
  con.mats <- plyr::alply(E, 3:4, consensus_matrix, .dims = TRUE) %>% 
    utils::relist(stats::setNames(replicate(
      dim(E)[4],
      list(structure(1:dim(E)[3], names = dimnames(E)[[3]]))),
      dimnames(E)[[4]]))
  con.cls <- mapply(function(cm, k) lapply(cm, function(x) hc(stats::dist(x), k = k)),
                    cm = con.mats, k = as.numeric(names(con.mats)),
                    SIMPLIFY = FALSE)
  out <- list(consensus_matrix = con.mats, consensus_class = con.cls) %>% 
    purrr::transpose() %>% 
    lapply(purrr::transpose)
  return(out)
}
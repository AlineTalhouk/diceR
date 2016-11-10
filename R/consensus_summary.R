#' Consensus summary
#' 
#' Given an object from \code{\link{ConClust}}, returns a list of consensus 
#' matrices and consensus classes for each clustering algorithm.
#' 
#' @param E output from \code{ConClust}
#' @return A list with summaries for each algorithm. Each algorithm has a list 
#'   with two elements: consensus_matrix and consensus_class
#' @family consensus functions
#' @author Derek Chiu
#' @export
#' @examples
#' data(hgsc)
#' dat <- t(hgsc[, -1])
#' x <- ConClust(dat, nc = 3:4, reps = 5, method = c("hcAEucl"))
#' cs <- consensus_summary(x)
#' str(cs)
consensus_summary <- function(E) {
  con.mats <- plyr::alply(E, 3:4, consensus_matrix, .dims = TRUE) %>% 
    utils::relist(stats::setNames(replicate(
      dim(E)[4],
      list(structure(1:dim(E)[3], names = dimnames(E)[[3]]))),
      dimnames(E)[[4]]))
  con.cls <- mapply(function(cm, k) lapply(cm, consensus_class, k = k),
                    cm = con.mats, k = as.numeric(names(con.mats)),
                    SIMPLIFY = FALSE)
  out <- list(consensus_matrix = con.mats, consensus_class = con.cls) %>% 
    purrr::transpose() %>% 
    lapply(purrr::transpose)
  return(out)
}
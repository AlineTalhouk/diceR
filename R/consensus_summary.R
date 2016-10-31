#' Consensus summary
#'
#' Given an object from \code{\link{ConClust}}, returns a list of consensus
#' matrices and consensus classes for each clustering algorithm.
#'
#' @param res output from \code{ConClust}
#' @param k number of clusters. The default is \code{NULL}, which returns 
#'   results for all k considered in \code{res}
#' @param progress logical; should a progress bar be displayed?
#' @return A list with summaries for each algorithm. Each algorithm has a list
#'   with two elements: consensus_matrix and consensus_class
#' @family consensus functions
#' @author Derek Chiu
#' @export
#' @examples
#' data(hgsc)
#' dat <- t(hgsc[, -1])
#' x <- ConClust(dat, nc = 2:4, reps = 5, method = c("hcAEucl"))
#' y1 <- consensus_summary(x)
#' y2 <- consensus_summary(x, k = 4)
#' str(y1)
#' str(y2)
consensus_summary <- function(res, k = NULL, progress = TRUE) {
  prog <- ifelse(progress, "text", "none")
  con.mats <- plyr::alply(res, 3:4, consensus_matrix, .progress = prog,
                          .dims = TRUE) %>% 
    utils::relist(stats::setNames(replicate(
      dim(res)[4],
      list(structure(1:dim(res)[3], names = dimnames(res)[[3]]))),
      dimnames(res)[[4]]))
  con.cls <- mapply(function(cm, k) lapply(cm, consensus_class, k = k),
                    cm = con.mats, k = as.numeric(names(con.mats)), SIMPLIFY = FALSE)
  out <- list(consensus_matrix = con.mats, consensus_class = con.cls) %>% 
    purrr::transpose() %>% 
    lapply(purrr::transpose)
  if (!is.null(k)) {
    out <- out[[as.character(k)]]
  }
  return(out)
}
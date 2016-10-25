#' Consensus summary
#'
#' Given an object from \code{\link{ConClust}}, returns a list of consensus
#' matrices and consensus classes for each clustering algorithm.
#'
#' @param res output from \code{ConClust}
#' @param progress logical; should a progress bar be displayed?
#' @param save logical; if \code{TRUE}, the returned object will be saved at
#'   each iteration as well as at the end.
#' @param file.name file name of the written object
#' @return A list with summaries for each algorithm. Each algorithm has a list
#'   with two elements: consensus_matrix and consensus_class
#' @family consensus functions
#' @author Derek Chiu
#' @importFrom utils relist
#' @export
#' @examples
#' data(hgsc)
#' dat <- t(hgsc[, -1])
#' x <- ConClust(dat, nc = 2:4, reps = 10, method = c("hcAEucl", "kmEucl"))
#' y1 <- consensus_summary(x)
#' y2 <- consensus_summary(x, k = 4)
#' str(y1)
#' str(y2)
consensus_summary <- function(res, k = NULL, progress = TRUE, save = FALSE,
                              file.name = "results_CC") {
  prog <- ifelse(progress, "text", "none")
  con.mats <- plyr::alply(res, 3:4, consensus_matrix, .progress = prog,
                          .dims = TRUE) %>% 
    relist(setNames(replicate(dim(res)[4],
                              list(structure(1:dim(res)[3],
                                             names = dimnames(res)[[3]]))),
                    dimnames(res)[[4]]))
  con.cls <- mapply(function(cm, k) lapply(cm, consensus_class, k = k),
                    cm = con.mats, k = as.numeric(names(con.mats)), SIMPLIFY = FALSE)
  out <- list(consensus_matrix = con.mats, consensus_class = con.cls) %>% 
    purrr::transpose() %>% 
    lapply(purrr::transpose)
  if (!is.null(k)) {
    out <- out[[as.character(k)]]
  }
  if (save) readr::write_rds(out, path = paste0(file.name, ".rds"),
                             compress = "xz")
  return(out)
}
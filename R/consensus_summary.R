#' Consensus summary
#'
#' Given an object from \code{\link{ConClust}}, returns a list of consensus
#' matrices and consensus classes for each clustering algorithm.
#'
#' @param res output from \code{ConClust}
#' @param k number of clusters to compute for consensus class
#' @param save logical; if \code{TRUE}, the returned object will be saved at
#'   each iteration as well as at the end.
#' @param file.name file name of the written object
#' @return A list with summaries for each algorithm. Each algorithm has a list
#'   with two elements: consensus_matrix and consensus_class
#' @family consensus functions
#' @author Derek Chiu
#' @export
#' @examples
#' data(hgsc)
#' dat <- t(hgsc[, -1])
#' x <- ConClust(dat, k = 4, reps = 10, method = "hcAEucl", save = FALSE)
#' y <- consensus_summary(x, k = 2)
#' str(y)
consensus_summary <- function(res, k, save = FALSE, file.name = "results_CC") {
  con.mats <- plyr::alply(res, 3, consensus_matrix, .progress = "text",
                          .dims = TRUE)
  con.cls <- plyr::llply(con.mats, consensus_class, k = k)
  z <- list(consensus_matrix = con.mats, consensus_class = con.cls)
  zv <- unlist(unname(z), recursive = FALSE)
  out <- split(setNames(zv, rep(names(z), lengths(z))), names(zv))
  if (save) readr::write_rds(out, path = paste0(file.name, ".rds"),
                             compress = "xz")
  return(out)
}

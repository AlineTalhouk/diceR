#' @param data data matrix with rows as samples and columns as variables 
#' @param cl.mat matrix of cluster assignments. Each row is an assignment for a
#'   different algorithm. Use \code{element = "class"} in
#'   \code{consensus_combine}.
#' @param cons.mat A list of consensus matrices, one for each algorithm. Use
#'   \code{element = "matrix"} in \code{consensus_combine}.
#' @param ref.cl reference class
#' @inheritParams consensus_combine
#' @return \code{consensus_compare} returns a data frame of the indices PAC and
#'   CHI in each column for each algorithm.
#' @rdname consensus_combine
#' @export
consensus_compare <- function(data, cl.mat, cons.mat, ref.cl = NULL,
                              alg.names = NULL) {
  if (!is.null(alg.names)) {
    an <- alg.names
  } else {
    an <- names(cons.mat)
  }
  x <- data.frame(data)
  ind.int <- data.frame(
    Algorithms = an,
    PAC = sapply(cons.mat, PAC, lower = 0.05, upper = 0.95),
    CHI = apply(cl.mat, 2, iv_chi, x = x),
    Compactness = apply(cl.mat, 2, iv_compactness, data = x)) %>% 
    cbind(plyr::aaply(cl.mat, 2, function(cl)
      unlist(iv_db_dunn(data = x, labels = cl)))) %>% 
    mutate_all(funs(structure(., names = an)))
  if (!is.null(ref.cl)) {
    ind.ext <- data.frame(
      Algorithms = an,
      Accuracy = apply(cl.mat, 2, ev_accuracy, ref.lab = ref.cl),
      plyr::aaply(cl.mat, 2, function(cl)
        unlist(ev_rand(pred.lab = cl, ref.lab = ref.cl))),
      NMI = apply(cl.mat, 2, ev_nmi, ref.lab = ref.cl)) %>% 
      mutate_all(funs(structure(., names = an)))
    return(list(internal = ind.int, external = ind.ext))
  } else {
    return(list(internal = ind.int))
  }
}
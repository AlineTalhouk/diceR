#' @param data data matrix with rows as samples and columns as variables
#' @param ref.cl reference class
#' @param plot logical; if \code{TRUE}, \code{graph_all} is called
#' @inheritParams consensus_combine
#' @return \code{consensus_evaluate} returns a data frame of the indices in each
#'   column for each algorithm.
#' @rdname consensus_combine
#' @export
consensus_evaluate <- function(data, k, ..., ref.cl = NULL, plot = TRUE) {
  x <- data.frame(data)
  cc.obj <- abind::abind(list(...), along = 3)
  cl.mat <- consensus_combine(cc.obj, k = k, element = "class")
  an <- dimnames(cc.obj)[3][[1]]
  ind.int <- data.frame(
    Algorithms = an,
    CHI = apply(cl.mat, 2, iv_chi, x = x),
    Compactness = apply(cl.mat, 2, iv_compactness, data = x)) %>% 
    cbind(plyr::aaply(cl.mat, 2, function(cl)
      unlist(iv_db_dunn(data = x, labels = cl))) %>% 
        rbind()) %>% 
    mutate_all(funs(structure(., names = an)))
  if (plot) {
    graph_all(cc.obj)
  }
  if (!is.null(ref.cl)) {
    ind.ext <- data.frame(
      Algorithms = an,
      Accuracy = apply(cl.mat, 2, ev_accuracy, ref.lab = ref.cl),
      plyr::aaply(cl.mat, 2, function(cl)
        unlist(ev_rand(pred.lab = cl, ref.lab = ref.cl))),
      NMI = apply(cl.mat, 2, ev_nmi, ref.lab = ref.cl)) %>% 
      mutate_all(funs(structure(., names = an)))
    return(list(k = k, internal = ind.int, external = ind.ext))
  } else {
    return(list(k = k, internal = ind.int))
  }
}
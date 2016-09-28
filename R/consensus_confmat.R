#' Confusion matrix between consensus clusters and true classes
#'
#' Compares the consensus class to the true class partitions.
#'
#' @param cl.cons vector of cluster assignments
#' @param cl.true true class labels
#' @param pred.lab label for the predicted clusters
#' @param ref.lab label for the reference classes
#' @return A confusion matrix with the predicted cluster assignments compared to
#'   the reference true class labels.
#' @family consensus functions
#' @author Derek Chiu
#' @export
#' @examples
#' set.seed(2)
#' a <- sample(1:4, 100, replace = TRUE)
#' b <- sample(1:4, 100, replace = TRUE)
#' consensus_confmat(a, b)
consensus_confmat <- function(cl.cons, cl.true, pred.lab = "Prediction",
                             ref.lab = "Reference") {
  pmat <- NULL
  confmat <- table(cl.cons, cl.true, dnn = c(eval(pred.lab), eval(ref.lab))) %>%
    min_fnorm() %>%
    use_series(pmat) %>%
    set_rownames(colnames(.)) %>%
    as.table()
  return(confmat)
}

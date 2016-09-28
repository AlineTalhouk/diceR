#' Classification Accuracy
#'
#' Given a sorted table, the accuracy is calculated based on proportion
#' correctly classified.
#'
#' The table must have rows and columns corresponding to the same classes. Then
#' the classification accuracy is computed as the sum of the diagonal entries
#' divided by the sum of all entries in the confusion matrix comparing the
#' predicted and reference classes.
#'
#' @param tbl table with predicted and reference classes correctly matched
#' @return classification accuracy, given as a proportion
#' @author Derek Chiu
#' @export
#' @examples
#' set.seed(1)
#' x <- matrix(rbinom(16, 20, 0.4), nrow = 4)
#' accuracy(x)
accuracy <- function(tbl) {
  return(sum(diag(tbl)) / sum(tbl))
}

#' Normalized Mutual Information
#'
#' Computes the NMI given a clustering assignment and true class labels.
#'
#' The function is adapted from the \code{mutinformation} function in the
#' \code{infotheo} package.
#' @param X clustering assignment
#' @param Y true class
#' @param method method of computing the entropy. Can be any one of "emp", "mm",
#'   "shrink", or "sg".
#' @return returns the normalized mutual information.
#' @author Derek Chiu
#' @references Strehl A, Ghosh J. Cluster ensembles: a knowledge reuse framework
#'   for combining multiple partitions. J. Mach. Learn. Res. 2002;3:583-617.
#' @export
#' @examples
#' set.seed(4)
#' X <- sample(1:4, 100, replace = TRUE)
#' Y <- sample(1:4, 100, replace = TRUE)
#' NMI(X, Y)
NMI <- function(X, Y, method = "emp") {
  U <- data.frame(Y, X)
  Hyx <- infotheo::entropy(U, method)
  Hx <- infotheo::entropy(X, method)
  Hy <- infotheo::entropy(Y, method)
  I <- ifelse(Hx + Hy - Hyx < 0, 0, Hx + Hy - Hyx)
  NMI <- I / sqrt(Hx * Hy)
  return(NMI)
}

#' Proportion of Ambiguous Clusters
#'
#' Given a consensus matrix, returns the proportion of ambiguous clusters (PAC).
#' This is a robust way to assess clustering performance.
#'
#' Since a consensus matrix is symmetric, we only look at its lower (or upper)
#' triangular matrix. The proportion of entries strictly between \code{lower}
#' and \code{upper} is the PAC. In a perfect clustering, the consensus matrix
#' would consist of only 0s and 1s, and the PAC assessed on the (0, 1) interval
#' would have a perfect score of 0. Using a (0.1, 0.9) interval for defining
#' ambiguity is common as well.
#'
#' @param cm consensus matrix. Should be symmetric and values between 0
#'   and 1.
#' @param lower the lower bound that determines what is ambiguous
#' @param upper the upper bound that determines what is ambiguous
#' @return the PAC is a score used in clustering performance. The lower it is
#'   the better, because we want minimal ambiguity amongst the consensus.
#' @author Derek Chiu
#' @export
#' @examples
#' set.seed(1)
#' x <- replicate(100, rbinom(100, 4, 0.2))
#' y <- consensus_matrix(x)
#' PAC(y, lower = 0.05, upper = 0.95)
PAC <- function(cm, lower = 0, upper = 1) {
  pac <- cm %>%
    extract(lower.tri(.)) %>%
    extract(and(is_greater_than(., lower), is_less_than(., upper))) %>%
    length() %>%
    divide_by(., length(cm[lower.tri(cm)]))
  return(pac)
}

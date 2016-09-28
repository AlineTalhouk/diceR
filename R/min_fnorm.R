#' Minimize Frobenius norm for between two matrices
#'
#' Finds a permutation of a matrix such that its Frobenius norm with another
#' matrix is minimized.
#'
#' Finds the permutation P of A such that \code{||PA - B||} is minimum in
#' Frobenius norm. Uses the linear-sum assignment problem (LSAP) solver in the
#' package \code{clue}. The default B is the identity matrix of same dimension,
#' so that the permutation of A maximizes its trace. This procedure is useful
#' for constructing a confusion matrix when we don't know the true class labels
#' of a predicted class and want to compare to a reference class.
#'
#' @param A data matrix we want to permute
#' @param B matrix whose distance with the permuted A we want to minimize. By
#'   default, \code{B <- diag(nrow(A))}, so the permutation maxmizes the trace
#'   of A.
#' @return Permuted matrix such that it is the permutation of A closest to B
#' @author Ravi Varadhan:
#'   https://stat.ethz.ch/pipermail/r-help/2010-April/236664.html
#' @export
#' @examples
#' A <- matrix(sample(1:25, size = 25, rep = FALSE), 5, 5)
#' min_fnorm(A)
min_fnorm <- function(A, B = diag(nrow(A))) {
  n <- nrow(A)
  D <- matrix(NA, n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      D[j, i] <- (sum((B[j, ] - A[i, ]) ^ 2))
    }
  }
  vec <- c(clue::solve_LSAP(D))
  return(list(pmat = A[vec, ], perm = vec, ord = order(vec)))
}

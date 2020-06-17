# magrittr placeholder
globalVariables(".")

#' Same binary operator as %in% except for partial matching
#' @noRd
`%pin%` <- function(x, table) pmatch(x, table, nomatch = 0L) > 0L

#' Minimize Frobenius norm for between two matrices
#'
#' Finds a permutation of a matrix such that its Frobenius norm with another
#' matrix is minimized.
#'
#' Finds the permutation P of A such that `||PA - B||` is minimum in Frobenius
#' norm. Uses the linear-sum assignment problem (LSAP) solver in the package
#' `clue`. The default B is the identity matrix of same dimension, so that the
#' permutation of A maximizes its trace. This procedure is useful for
#' constructing a confusion matrix when we don't know the true class labels of a
#' predicted class and want to compare to a reference class.
#'
#' @param A data matrix we want to permute
#' @param B matrix whose distance with the permuted A we want to minimize. By
#'   default, `B <- diag(nrow(A))`, so the permutation maximizes the trace of A.
#' @return Permuted matrix such that it is the permutation of A closest to B
#' @author Ravi Varadhan:
#'   https://stat.ethz.ch/pipermail/r-help/2010-April/236664.html
#' @export
#' @examples
#'
#' set.seed(1)
#' A <- matrix(sample(1:25, size = 25, rep = FALSE), 5, 5)
#' min_fnorm(A)
min_fnorm <- function(A, B = diag(nrow(A))) {
  n <- nrow(A)
  D <- matrix(NA, n, n)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      D[j, i] <- sum((B[j, ] - A[i, ]) ^ 2)
    }
  }
  vec <- clue::solve_LSAP(D)
  list(pmat = A[vec, ], perm = vec, ord = order(vec))
}

#' Relabel classes to a standard
#'
#' Relabel clustering categories to match to a standard by minimizing the
#' Frobenius norm between the two labels.
#'
#' @param pred.cl vector of predicted cluster assignments
#' @param ref.cl vector of reference labels to match to
#' @return A vector of relabeled cluster assignments
#' @author Aline Talhouk
#' @export
#' @examples
#' set.seed(2)
#' pred <- sample(1:4, 100, replace = TRUE)
#' true <- sample(1:4, 100, replace = TRUE)
#' relabel_class(pred, true)
relabel_class <- function(pred.cl, ref.cl) {
  pred.cl %>%
    factor(levels = sort(unique(ref.cl))) %>%
    table(., ref.cl) %>%
    min_fnorm() %>%
    magrittr::use_series("perm") %>%
    factor(pred.cl, levels = ., labels = levels(factor(ref.cl))) %>%
    as.integer()
}

#' Helper function for k_modes and majority_voting to flatten (and relabel) E
#' @noRd
flatten_E <- function(E, is.relabelled) {
  # take E imputed and reshape into a flat matrix
  if (length(dim(E)) > 2) {
    flat_E <- matrix(E, nrow = dim(E)[1], ncol = prod(dim(E)[2:3]))
  } else {
    flat_E <- E
  }
  # relabel using first clustering as reference
  if (!is.relabelled & ncol(flat_E) > 1) {
    flat_E <- cbind(flat_E[, 1],
                    apply(flat_E[, -1, drop = FALSE], 2, relabel_class,
                          ref.cl = flat_E[, 1]))
  }
  flat_E
}

#' Column min for matrix
#' @param x a matrix
#' @param na.rm logical; Should missing values be omitted from consideration?
#' @noRd
colMin <- function(x, na.rm = TRUE) {
  assertthat::assert_that(is.matrix(x), is.numeric(x))
  apply(x, 2, min, na.rm = na.rm)
}

#' Flattened row indices for upper triangular matrix elements of diag(n)
#' @noRd
upper_tri_row <- function(n) unlist(lapply(seq_len(n - 1), seq_len))

#' Flattened column indices for upper triangular matrix elements of diag(n)
#' @noRd
upper_tri_col <- function(n) rep(seq_len(n)[-1], seq_len(n - 1))

#' Which rows of matrix x contain the value?
#' @noRd
which_row <- function(x, value) which(x == value, arr.ind = TRUE)[, "row"]

#' Check if a single number is a positive integer
#' @noRd
is_pos_int <- function(x) {
  assertthat::assert_that(length(x) == 1)
  if (x <= 0) {
    return(FALSE)
  } else {
    if (x %% 1 == 0) {
      return(TRUE)
    } else{
      return(FALSE)
    }
  }
}

#' Redirect any console printouts from print() or cat() to null device
#' @references
#'   http://stackoverflow.com/questions/5310393/disabling-the-cat-command
#' @noRd
sink_output <- function(expr) {
  tmp <- tempfile()
  sink(tmp)
  on.exit(sink())
  on.exit(file.remove(tmp), add = TRUE)
  invisible(force(expr))
}

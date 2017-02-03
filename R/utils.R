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
#' Finds the permutation P of A such that \code{||PA - B||} is minimum in
#' Frobenius norm. Uses the linear-sum assignment problem (LSAP) solver in the
#' package \code{clue}. The default B is the identity matrix of same dimension,
#' so that the permutation of A maximizes its trace. This procedure is useful
#' for constructing a confusion matrix when we don't know the true class labels
#' of a predicted class and want to compare to a reference class.
#' 
#' @param A data matrix we want to permute
#' @param B matrix whose distance with the permuted A we want to minimize. By
#' default, \code{B <- diag(nrow(A))}, so the permutation maxmizes the trace of
#' A.
#' @return Permuted matrix such that it is the permutation of A closest to B
#' @author Ravi Varadhan:
#' https://stat.ethz.ch/pipermail/r-help/2010-April/236664.html
#' @export
#' @examples
#' 
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

#' Relabel classes to a standard
#'
#' Relabel clustering categories to match to a standard by minimizing
#' the Frobenius norm between the two labels.
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
  perm <- pred.cl %>% 
    factor(levels = sort(unique(ref.cl))) %>% 
    table(., ref.cl) %>%
    min_fnorm() %>%
    use_series(perm)
  res <- factor(pred.cl, levels = perm, labels = levels(factor(ref.cl))) %>% 
    as.integer()
  return(res)
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
                    apply(flat_E[, -1], 2, relabel_class, ref.cl = flat_E[, 1]))
  }
  return(flat_E)
}

#' Column min for matrix
#' @param x a matrix
#' @param na.rm logical; Should missing values be omitted from consideration?
#' @noRd
colMin <- function(x, na.rm = TRUE) {
  assertthat::assert_that(is.matrix(x), is.numeric(x))
  return(apply(x, 2, min, na.rm = na.rm))
}

#' Find coordinates of matrix x that give element n
#' @noRd
coord <- function(x, n) {
  assertthat::assert_that(is.matrix(x), is.numeric(x), is.numeric(n),
                          n %in% x, length(n) == 1)
  res <- which(x == n, arr.ind = TRUE)
  return(stats::setNames(unlist(apply(res, 2, list), recursive = FALSE),
                         c("rows", "cols")))
}

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

#' function to sort rows of a matrix
#' @noRd
sortMatrixRowWise <- function(M, order = c("ascending", "descending")) {
  assertthat::assert_that(nrow(M) >= 1, is.matrix(M), is.numeric(M))
  for (i in 1:nrow(M)) {
    temp <- M[i, ]
    temp <- switch(match.arg(order),
                   ascending = sort(temp),
                   descending = sort(temp, decreasing = TRUE))
    M[i, ] <- temp
  }
  return(M)
}
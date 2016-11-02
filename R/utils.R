# magrittr placeholder
globalVariables(".")

#' Custom hook functions for DIvisive ANAlysis clustering algorithm
#'
#' @param d distance matrix
#' @param k scalar indicating number of clusters to cut tree into
#' @noRd
diana_hook <- function(d, k) {
  tmp <- cluster::diana(d, diss = TRUE)
  a <- stats::cutree(tmp, k)
  return(a)
}

#' Prepare data for consensus clustering
#'
#' Remove variables with low signal and scale before consensus clustering
#'
#' The \code{min.sd} argument is used to filter the feature space for only
#' highly variable features. Only features with a standard deviation across all
#' samples greater than \code{min.sd} will be used.
#'
#' @param data data matrix with rows as samples and columns as variables
#' @param min.sd minimum standard deviation threshold. See details.
#' @return dataset prepared for usage in \code{ConClust}
#' @author Derek Chiu
#' @export
#' @examples
#' set.seed(2)
#' x <- replicate(10, rnorm(100))
#' prepare_data(x)
prepare_data <- function(data, min.sd = 1) {
  dat.out <- data %>%
    magrittr::extract(apply(., 1, function(x) !any(is.na(x))),
                      apply(., 2, function(x) stats::sd(x, na.rm = TRUE)) >
                        min.sd) %>%
    scale()
  return(dat.out)
}

#' Calculate pairwise Spearman correlational distances using
#' bioDist::spearman.dist defaults
#' @references https://github.com/Bioconductor-mirror/bioDist/blob/master/R/spearman.dist.R
#' @noRd
spearman_dist <- function(x) {
  rvec <- cor(t(x), method = "spearman") %>% 
    abs() %>% 
    magrittr::subtract(1, .) %>% 
    magrittr::extract(lower.tri(.))
  attributes(rvec) <- list(Size = nrow(x), Labels = rownames(x), Diag = FALSE,
                           Upper = FALSE, methods = "spearman", class = "dist")
  return(rvec)
}

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
#' @examples
#' 
#' A <- matrix(sample(1:25, size = 25, rep = FALSE), 5, 5)
#' min_fnorm(A)
#' 
#' @export min_fnorm
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
#' @param cl.pred vector of predicted cluster assignments
#' @param cl.ref vector of reference labels to match to
#' @return A vector of relabeled cluster assignments
#' @author Aline Talhouk
#' @export
#' @examples
#' set.seed(2)
#' pred <- sample(1:4, 100, replace = TRUE)
#' true <- sample(1:4, 100, replace = TRUE)
#' relabel_class(pred, true)
relabel_class <- function(cl.pred, cl.ref) {
  perm <- table(cl.pred, cl.ref) %>%
    min_fnorm() %>%
    use_series(perm)
  res <- factor(cl.pred, levels = perm, labels = levels(factor(cl.ref)))
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
  if (!is.relabelled) {
    flat_E <- cbind(flat_E[, 1],
                    apply(flat_E[, -1], 2, relabel_class, cl.ref = flat_E[, 1]))
  }
  return(flat_E)
}

#' Form Row and Column Extremes
#' 
#' Form row and column maxs and mins for matrices
#' 
#' @param x a matrix
#' @param na.rm logical; Should missing values be omitted from consideration?
#'   
#' @return a vector containing the maximums or minimums for each column in
#'   \code{x}
#' @export
#' 
#' @examples
#' a <- matrix(1:9, ncol = 3)
#' colMax(a)
#' colMin(a)
colMax <- function(x, na.rm = TRUE) {
  assertthat::assert_that(is.matrix(x), is.numeric(x))
  return(apply(x, 2, max, na.rm = na.rm))
}

#' @rdname colMax
#' @export
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
sortMatrixRowWise <- function(M, order) {
  assertthat::assert_that(nrow(M) >= 1, is.matrix(M), is.numeric(M))
  for (i in 1:nrow(M)) {
    temp <- M[i, ]
    if (order == "ascending") {
      temp <- sort(temp)
    } else if (order == "descending") {
      temp <- sort(temp, decreasing = TRUE)
    } else {
      stop("Invalid order argument in sortMatrixRowWise")
    }
    M[i, ] <- temp
  }
  return(M)
}
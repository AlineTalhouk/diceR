#' Create hierarchical cluster tree
#'
#' @param Y a distance vector
#' @param method linkage method; can be "average", "single", or "complete"
#' @references function linkage in MATLAB for average, single, and complete linkage
#' @return matrix denoting hierarchical cluster tree
#' @export
#'
#' @examples
#' data("E_LCE")
#' linkage(stod(asrs(E_LCE, 0.8)), "average")
linkage <- function(Y, method) {
  assertthat::assert_that(is.vector(Y))
  assertthat::assert_that(is.numeric(Y))
  assertthat::assert_that(method %in% c("complete", "average", "single"))
  n <- length(Y)
  m <- ceiling(sqrt(2 * n))
  Z <- pracma::zeros(m - 1, 3)
  N <- rep(0, 2 * m - 1)
  N[1:m] <- 1
  n <- m
  R <- 1:n
  for (s in 1:(n - 1)) {
    if (method == "average") {
      if ((m - 1) >= 2) {
        p <- seq(m - 1, 2, -1)
      }
      I <- rep(0, 0.5 * m * (m - 1))
      I[cumsum(c(1, p))] <- 1
      I <- cumsum(I)
      J <- rep(1, 0.5 * m * (m - 1))
      J[cumsum(p) + 1] <- 2 - p
      J[1] <- 2
      J <- cumsum(J)
      W <- N[R[I]] * N[R[J]]
      W <- W[!is.na(W)]
      tempVec <- Y / W
      v <- min(tempVec)
      k <- which(tempVec == min(tempVec))[1]
    } else{
      v <- min(Y)
      k <- which(Y == min(Y))[1]
    }
    if ((length(m) != 0) & (length(k) != 0)) {
      i <- floor(m + 1 / 2 - sqrt(m ^ 2 - m + 1 / 4 - 2 * (k - 1)))
      j <- k - (i - 1) * (m - i / 2) + i
    }
    Z[s, ] <- cbind(R[i], R[j], v)
    tempUIJ <- UIJ(i, j, m)
    U <- tempUIJ$U
    I <- tempUIJ$I
    J <- tempUIJ$J
    
    if (!is.null(I)) {
      if (method == "single") {
        Y[I] <- colMin(rbind(Y[I], Y[J]))
      } else if (method == "complete") {
        Y[I] <- colMax(rbind(Y[I], Y[J]))
      } else if (method == "average") {
        Y[I] <- Y[I] + Y[J]
      } else{
        stop("Invalid method for function linkage")
      }
    }
    J <- c(J, i * (m - (i + 1) / 2) - m + j)
    Y <- Y[-J]
    m <- m - 1
    N[n + s] <- N[R[i]] + N[R[j]]
    R[i] <- n + s
    if (((n - 1) >= j) & (n >= (j + 1))) {
      R[j:(n - 1)] <- R[(j + 1):n]
    }
  }
  Z[, 1:2] <- sortMatrixRowWise(Z[, 1:2], "ascending")
  return(Z)
}

#' Update U, I, and J for \code{linkage}
#'
#' @param i see i in \code{linkage}
#' @param j see j in \code{linkage}
#' @param m see m in \code{linkage}
#' @noRd
UIJ <- function(i, j, m) {
  assertthat::assert_that(is.numeric(i))
  assertthat::assert_that(is.numeric(j))
  assertthat::assert_that(is.numeric(m))
  assertthat::assert_that(length(i) == 1)
  assertthat::assert_that(length(j) == 1)
  assertthat::assert_that(length(m) == 1)
  I1 <- NULL
  I2 <- NULL
  I3 <- NULL
  I <- NULL
  J <- NULL
  if (i >= 2) {
    I1 <- 1:(i - 1)
  }
  if ((j - 1) >= (i + 1)) {
    I2 <- (i + 1):(j - 1)
  }
  if (m >= (j + 1)) {
    I3 <- (j + 1):m
  }
  if (!is.null(I1)) {
    I <- c(I, I1 * (m - (I1 + 1) / 2) - m + i)
    J <- c(J, I1 * (m - (I1 + 1) / 2) - m + j)
  }
  if (!is.null(I2)) {
    I <- c(I, i * (m - (i + 1) / 2) - m + I2)
    J <- c(J, I2 * (m - (I2 + 1) / 2) - m + j)
  }
  if (!is.null(I3)) {
    I <- c(I, i * (m - (i + 1) / 2) - m + I3)
    J <- c(J, j * (m - (j + 1) / 2) - m + I3)
  }
  U <- c(I1, I2, I3)
  return(list(U = U, I = I, J = J))
}
#' Similarity Matrices
#' 
#' \code{srs} computes the simrank based similarity matrix, \code{asrs} computes
#' the approximated simrank based similarity matrix, and \code{cts} computes the
#' connected triple based similarity matrix
#' 
#' @param E an N by M matrix of cluster ensembles
#' @param dc decay factor, ranges from 0 to 1 inclusive
#' @param R number of iterations for simrank algorithm
#' @references MATLAB function srs, cts, asrs in package LinkCluE by Simon Garrett   
#' @return an N by N SRS, ASRS, or CTS matrix
#' @name similarity
#' @author Johnson Liu
#' @export
#' 
#' @examples 
#' set.seed(1)
#' E <- matrix(rep(sample(1:4, 800, replace = TRUE)), nrow = 100)
#' SRS <- srs(E = E, dc = 0.8, R = 3)
#' ASRS <- asrs(E = E, dc = 0.8)
#' CTS <- cts(E = E, dc = 0.8)
srs <- function(E, dc, R) {
  assertthat::assert_that(is.matrix(E), is.numeric(R), dc >= 0 && dc <= 1,
                          is_pos_int(R) == TRUE)
  n <- nrow(E)
  M <- ncol(E)
  E.new <- relabel_clusters(E)
  E <- E.new$newE
  no_allcl <- E.new$no_allcl
  S <- diag(x = 1, nrow = n, ncol = n)
  C <- diag(x = 1, nrow = no_allcl, ncol = no_allcl)
  
  for (r in 1:(R - 1)) {
    S1 <- diag(x = 1, nrow = n, ncol = n)
    for (i in 1:(n - 1)) {
      Ni <- E[i, ]
      for (ii in (i + 1):n) {
        sum_sim <- 0
        Nii <- E[ii, ]
        for (k in 1:M) {
          for (kk in 1:M) {
            sum_sim <- sum_sim + C[Ni[k], Nii[kk]]
          }
        }
        S1[i, ii] <- (dc / (M * M)) * sum_sim
      }
    }
    S1 <- S1 + t(S1)
    for (i in 1:n) {
      S1[i, i] <- 1
    }
    C1 <- diag(x = 1, nrow = no_allcl, ncol = no_allcl)
    for (i in 1:(no_allcl - 1)) {
      Ni <- coord(E, i)$rows
      col <- coord(E, i)$cols
      nki <- length(Ni)
      for (ii in (i + 1):no_allcl) {
        sum_sim <- 0
        Nii <- coord(E, ii)$rows
        col <- coord(E, ii)$cols
        nkii <- length(Nii)
        for (k in 1:nki) {
          for (kk in 1:nkii) {
            sum_sim <- sum_sim + S[Ni[k], Nii[kk]]
          }
        }
        if ((nki * nkii) > 0) {
          C1[i, ii] <- dc * sum_sim / (nki * nkii)
        }
      }
    }
    C1 <- C1 + t(C1)
    for (i in 1:no_allcl) {
      C1[i, i] <- 1
    }
    S <- S1
    C <- C1
  }
  return(S)
}

#' @rdname similarity
#' @export
asrs <- function(E, dc) {
  assertthat::assert_that(is.matrix(E), is.numeric(E), dc >= 0 && dc <= 1)
  n <- nrow(E)
  M <- ncol(E)
  E.new <- relabel_clusters(E)
  E <- E.new$newE
  no_allcl <- E.new$no_allcl
  wcl <- weigh_clusters(E)
  CS <- matrix(rep(0, no_allcl * no_allcl), nrow = no_allcl)
  for (i in 1:(no_allcl - 1)) {
    Ni <- wcl[i, ]
    ni <- length(Ni[which(Ni > 0)])
    for (j in (i + 1):no_allcl) {
      Nj <- wcl[j, ]
      nj <- length(Nj[which(Nj > 0)])
      if ((ni * nj) > 0) {
        CS[i, j] <- (Ni %*% Nj) / (ni * nj)
      }
    }
  }
  if (max(CS) > 0) {
    CS <- CS / max(CS)
  }
  CS <- CS + t(CS)
  for (i in 1:no_allcl) {
    CS[i, i] <- 1
  }
  S <- matrix(rep(0, n * n), nrow = n)
  for (i in 1:(n - 1)) {
    for (ii in (i + 1):n) {
      for (j in 1:M) {
        for (jj in 1:M) {
          if (CS[E[i, j], E[ii, jj]] == 1) {
            S[i, ii] <- S[i, ii] + 1
          } else{
            S[i, ii] <- S[i, ii] + dc * CS[E[i, j], E[ii, jj]]
          }
        }
      }
    }
  }
  S <- S / (M * M)
  S <- S + t(S)
  for (i in 1:n) {
    S[i, i] <- 1
  }
  return(S)
}

#' @rdname similarity
#' @export
cts <- function(E, dc) {
  assertthat::assert_that(is.matrix(E), dc >= 0 && dc <= 1)
  n <- nrow(E)
  M <- ncol(E)
  E.new <- relabel_clusters(E)
  E <- E.new$newE
  no_allcl <- E.new$no_allcl
  wcl <- weigh_clusters(E)
  wCT <- matrix(rep(0, no_allcl * no_allcl), nrow = no_allcl)
  maxCl <- apply(E, 2, max)
  minCl <- apply(E, 2, min)
  for (q in 1:M) {
    for (i in minCl[q]:(maxCl[q] - 1)) {
      Ni <- wcl[i, ]
      for (j in (i + 1):(maxCl[q])) {
        Nj <- wcl[j, ]
        wCT[i, j] <- sum(colMin(rbind(Ni, Nj)))
      }
    }
  }
  if (max(wCT) > 0) {
    wCT <- wCT / max(wCT)
  }
  wCT <- wCT + t(wCT)
  for (i in 1:no_allcl) {
    wCT[i, i] <- 1
  }
  S <- matrix(rep(0, n * n), nrow = n)
  for (m in 1:M) {
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        if (E[i, m] == E[j, m]) {
          S[i, j] = S[i, j] + 1
        } else{
          S[i, j] = S[i, j] + dc * wCT[E[i, m], E[j, m]]
        }
      }
    }
  }
  S <- S / M
  S <- S + t(S)
  for (i in 1:n) {
    S[i, i] <- 1
  }
  return(S)
}
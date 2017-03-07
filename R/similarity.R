#' Similarity Matrices
#' 
#' \code{srs} computes the simrank based similarity matrix, \code{asrs} computes
#' the approximated simrank based similarity matrix, and \code{cts} computes the
#' connected triple based similarity matrix.
#' 
#' @param E an N by M matrix of cluster ensembles
#' @param dc decay factor, ranges from 0 to 1 inclusive
#' @param R number of iterations for \code{srs}
#' @return an N by N SRS, ASRS, or CTS matrix
#' @name similarity
#' @author Johnson Liu
#' @references MATLAB functions srs, cts, asrs in package LinkCluE by Simon
#'   Garrett
#' @export
#' 
#' @examples 
#' set.seed(1)
#' E <- matrix(rep(sample(1:4, 800, replace = TRUE)), nrow = 100)
#' SRS <- srs(E = E, dc = 0.8, R = 3)
#' ASRS <- asrs(E = E, dc = 0.8)
#' CTS <- cts(E = E, dc = 0.8)
#' str_out <- sapply(list(SRS, ASRS, CTS), str)
srs <- function(E, dc, R) {
  assertthat::assert_that(is.matrix(E), is.numeric(E),
                          is.numeric(R), is_pos_int(R), dc >= 0 && dc <= 1)
  n <- nrow(E)
  M <- ncol(E)
  E.new <- relabel_clusters(E)
  E <- E.new$newE
  no_allcl <- E.new$no_allcl
  S <- diag(x = 1, nrow = n, ncol = n)
  C <- diag(x = 1, nrow = no_allcl, ncol = no_allcl)
  
  for (r in seq_len(R - 1)) {
    S1 <- diag(x = 1, nrow = n, ncol = n) %>% 
      inset(upper.tri(.), purrr::map2_dbl(
        upper_tri_row(n), upper_tri_col(n),
        ~ (dc / (M * M)) * sum(C[E[.x, ], E[.y, ]]))) %>% 
      add(t(.)) %>% 
      inset(row(.) == col(.), 1)
    C1 <- diag(x = 1, nrow = no_allcl, ncol = no_allcl) %>% 
      inset(upper.tri(.), purrr::map2_dbl(
        upper_tri_row(no_allcl), upper_tri_col(no_allcl), ~ {
          Ni <- which_row(E, .x)
          Nii <- which_row(E, .y)
          if (length(Ni) * length(Nii)) {
            dc * sum(S[Ni, Nii]) / (length(Ni) * length(Nii))
          }
        }
      )) %>% 
      add(t(.)) %>% 
      inset(row(.) == col(.), 1)
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
  CS <- sapply(seq_len(no_allcl), function(j) {
    sapply(seq_len(no_allcl), function(i) {
      if (i < no_allcl & i < j) {
        Ni <- wcl[i, ]
        ni <- length(Ni[which(Ni > 0)])
        Nj <- wcl[j, ]
        nj <- length(Nj[which(Nj > 0)])
        if ((ni * nj > 0))
          return((Ni %*% Nj) / (ni * nj))
      } else {
        return(0)
      }
    })
  })
  if (max(CS) > 0) {
    CS <- CS / max(CS)
  }
  CS <- CS + t(CS)
  CS[row(CS) == col(CS)] <- 1
  S <- matrix(rep(0, n * n), nrow = n)
  for (i in 1:(n - 1)) {
    for (ii in (i + 1):n) {
      for (j in 1:M) {
        for (jj in 1:M) {
          if (CS[E[i, j], E[ii, jj]] == 1) {
            S[i, ii] <- S[i, ii] + 1
          } else {
            S[i, ii] <- S[i, ii] + dc * CS[E[i, j], E[ii, jj]]
          }
        }
      }
    }
  }
  S <- S / (M * M)
  S <- S + t(S)
  S[row(S) == col(S)] <- 1
  return(S)
}

#' @rdname similarity
#' @export
cts <- function(E, dc) {
  assertthat::assert_that(is.matrix(E), is.numeric(E), dc >= 0 && dc <= 1)
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
  wCT[row(wCT) == col(wCT)] <- 1
  S <- matrix(rep(0, n * n), nrow = n)
  for (m in 1:M) {
    for (i in 1:(n - 1)) {
      for (j in (i + 1):n) {
        if (E[i, m] == E[j, m]) {
          S[i, j] = S[i, j] + 1
        } else {
          S[i, j] = S[i, j] + dc * wCT[E[i, m], E[j, m]]
        }
      }
    }
  }
  S <- S / M
  S <- S + t(S)
  S[row(S) == col(S)] <- 1
  return(S)
}

#' Relabel clusters in ensemble \code{E}
#'
#' @param E N by M cluster ensemble matrix
#' @return A list with elements
#' \item{newE}{N by M relabelled cluster ensemble matrix}
#' \item{no_allcl}{total number of clusters in the ensemble}
#' @author Johnson Liu
#' @references MATLAB function relabelCl by Simon Garrett in LinkCluE package
#' @noRd
relabel_clusters <- function(E) {
  if (!all(apply(E, 1:length(dim(E)), is_pos_int)))
    stop("Error: one of the entries in the input matrix is not a positive integer.")
  N <- nrow(E)
  M <- ncol(E)
  newE <- matrix(rep(0, N * M), nrow = N)
  ucl <- sort(unique(E[, 1]))
  if (max(E[, 1] != length(ucl)) == 1) {
    for (j in 1:length(ucl)) {
      newE[which(E[, 1] == ucl[j]), 1] <- j
    }
  }
  for (i in 2:M) {
    ucl <- sort(unique(E[, i]))
    prevCl <- length(sort(unique(c(newE[, 1:(i - 1)]))))
    for (j in 1:length(ucl)) {
      newE[E[, i] == ucl[j], i] <- prevCl + j
    }
  }
  no_allcl <- max(newE)
  return(list(no_allcl = no_allcl, newE = newE))
}

#' Compute weight for each pair of clusters using their shared members (Jaccard
#' coefficient)
#' 
#' @param E N by M cluster ensemble matrix
#' @return a p by p weighted cluster matrix where p denotes number of classes
#' @author Johnson Liu
#' @references MATLAB function weightCl by Simon Garrett in package LinkCluE   
#' @noRd
weigh_clusters <- function(E) {
  if (!all(apply(E, 1:length(dim(E)), is_pos_int)))
    stop("Error: one of the entries in the input matrix is not a positive integer.")
  N <- nrow(E)
  no_allcl <- max(E)
  pc <- matrix(rep(0, N * no_allcl), nrow = N)
  for (i in 1:N) {
    pc[i, E[i, ]] <- 1
  }
  wcl <- matrix(rep(0, no_allcl ^ 2), nrow = no_allcl)
  for (i in 1:(no_allcl - 1)) {
    for (j in (i + 1):no_allcl) {
      tmp <- sum(as.numeric((pc[, i] + pc[, j])) > 0)
      if (tmp > 0) {
        wcl[i, j] <- sum(as.numeric((pc[, i] + pc[, j])) == 2) / tmp
      }
    }
  }
  wcl <- wcl + t(wcl)
  return(wcl)
}
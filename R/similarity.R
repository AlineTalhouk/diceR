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
#' @author Johnson Liu, Derek Chiu
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
  S <- diag(1, n)
  C <- diag(1, no_allcl)
  
  for (r in seq_len(R - 1)) {
    S1 <- diag(1, n) %>% 
      inset(upper.tri(.), purrr::map2_dbl(
        upper_tri_row(n), upper_tri_col(n),
        ~ (dc / (M * M)) * sum(C[E[.x, ], E[.y, ]]))) %>% 
      add(t(.)) %>% 
      inset(row(.) == col(.), 1)
    C1 <- diag(1, no_allcl) %>% 
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
  CS <- vapply(seq_len(no_allcl), function(j) {
    vapply(seq_len(no_allcl), function(i) {
      if (i < no_allcl & i < j) {
        Ni <- wcl[i, ]
        ni <- length(Ni[Ni > 0])
        Nj <- wcl[j, ]
        nj <- length(Nj[Nj > 0])
        if (ni * nj)
          return((Ni %*% Nj) / (ni * nj))
      } else {
        return(0)
      }
    }, double(1))
  }, double(no_allcl)) %>% 
    inset(max(.) > 0, . / max(.)) %>% 
    add(t(.)) %>% 
    inset(row(.) == col(.), 1)
  S <- diag(0, n) %>% 
    inset(upper.tri(.), purrr::map(seq_len(n)[-1], function(i) {
      purrr::map_dbl(seq_len(i - 1), function(ii) {
        cse <- CS[E[i, ], E[ii, ]]
        sum(dc * cse[cse != 1]) + sum(cse == 1)
      })
    }) %>% 
      purrr::flatten_dbl()) %>% 
    divide_by(M * M) %>% 
    add(t(.)) %>% 
    inset(row(.) == col(.), 1)
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
  ind <- seq_len(no_allcl) %>% 
    split(ceiling(. / (no_allcl / M))) %>% 
    purrr::map(utils::combn, 2) %>% 
    do.call(cbind, .) %>% 
    t()
  wCT <- diag(0, no_allcl) %>% 
    inset(ind, ind %>% 
            purrr::array_branch(margin = 2) %>% 
            purrr::pmap_dbl(~ sum(colMin(wcl[c(.x, .y), ])))) %>% 
    inset(max(.) > 0, . / max(.)) %>% 
    add(t(.)) %>% 
    inset(row(.) == col(.), 1)
  S <- diag(0, n) %>% 
    inset(upper.tri(.), purrr::map(seq_len(n)[-1], function(i) {
      purrr::map_dbl(seq_len(i - 1), function(j) {
        Ei <- E[i, ]
        Ej <- E[j, ]
        sum(dc * wCT[Ei, Ej][, Ei != Ej]) + sum(wCT[Ei, Ej][, Ei == Ej])
      })
    }) %>% 
      purrr::flatten_dbl()) %>% 
    divide_by(M) %>% 
    add(t(.)) %>% 
    inset(row(.) == col(.), 1)
  return(S)
}

#' Relabel clusters in ensemble \code{E}
#'
#' @param E N by M cluster ensemble matrix
#' @return A list with elements
#' \item{newE}{N by M relabelled cluster ensemble matrix}
#' \item{no_allcl}{total number of clusters in the ensemble}
#' @author Johnson Liu, Derek Chiu
#' @references MATLAB function relabelCl by Simon Garrett in LinkCluE package
#' @noRd
relabel_clusters <- function(E) {
  if (!all(vapply(as.vector(E), is_pos_int, logical(1))))
    stop("Error: one of the entries in the input matrix is not a positive integer.")
  N <- nrow(E)
  M <- ncol(E)
  newE <- matrix(0, nrow = N, ncol = M)
  for (i in seq_len(M)) {
    ucl <- sort(unique(E[, i]))
    prevCl <- n_distinct(c(newE[, seq_len(i - 1)]))
    for (j in seq_along(ucl)) {
      newE[E[, i] == ucl[j], i] <- prevCl + j
    }
  }
  return(list(no_allcl = max(newE), newE = newE))
}

#' Compute weight for each pair of clusters using their shared members (Jaccard
#' coefficient)
#' 
#' @param E N by M cluster ensemble matrix
#' @return a p by p weighted cluster matrix where p denotes number of classes
#' @author Johnson Liu, Derek Chiu
#' @references MATLAB function weightCl by Simon Garrett in package LinkCluE   
#' @noRd
weigh_clusters <- function(E) {
  if (!all(vapply(as.vector(E), is_pos_int, logical(1))))
    stop("Error: one of the entries in the input matrix is not a positive integer.")
  N <- nrow(E)
  no_allcl <- max(E)
  pc <- matrix(0, nrow = N, ncol = no_allcl)
  for (i in seq_len(N)) {
    pc[i, E[i, ]] <- 1
  }
  wcl <- diag(0, no_allcl) %>% 
    inset(upper.tri(.), purrr::map2_dbl(
      upper_tri_row(no_allcl), upper_tri_col(no_allcl), ~ {
        tmp <- pc[, .x] + pc[, .y]
        ifelse(sum(tmp) > 0, sum(tmp == 2) / sum(tmp > 0), 0)
      }
    )) %>% 
    add(t(.))
  return(wcl)
}
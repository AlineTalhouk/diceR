#' @param quantile if an algorithm's sum of ranks across indices is in the 
#'   \code{quantile} quantile, the algorithm is kept. Otherwise it is removed 
#'   (trimmed).
#' @param reweigh logical; if \code{TRUE}, after trimming out poor performing
#'   algorithms, each algorithm is reweighed depending on its internal indices.
#' @param show.eval logical; if \code{TRUE} (default), show the evaluation 
#'   output from \code{consensus_evaluate}
#' @inheritParams consensus_evaluate
#' @return \code{consensus_trim} returns a list with three elements 
#'   \item{alg.keep}{algorithms kept}
#'   \item{alg.remove}{algorithms removed}
#'   \item{eval}{if \code{show.eval = TRUE}, the evaluation output is returned, 
#'   otherwise \code{NULL}}
#'   \item{data.new}{A new version of a \code{consensus_cluster} data object.
#'   Potentially no trimming depending on \code{quantile}} value.
#' @rdname consensus_combine
#' @export
consensus_trim <- function(data, ..., cons.cl = NULL, ref.cl = NULL,
                           quantile = 0.75, reweigh = FALSE, show.eval = TRUE) {
  # Evaluate and extract chosen k
  Sumrank <- Quantile <- NULL
  cc.obj <- abind::abind(list(...), along = 3)
  z <- consensus_evaluate(data = data, cc.obj, cons.cl = cons.cl,
                          ref.cl = ref.cl, plot = FALSE)
  k <- as.character(min(z$k))
  zk <- z$internal[[k]]
  alg.all <- zk$Algorithms
  
  # Separate algorithms into those from clusterCrit (main), and (others)
  z.main <- zk %>% 
    extract(!names(.) %in% c("Algorithms", "Compactness", "Connectivity") &
              purrr::map_lgl(., ~ all(!is.nan(.x))))
  z.other <- zk %>% 
    extract(c("Compactness", "Connectivity"))
  
  # Which algorithm is the best for each index?
  bests <- purrr::map2_int(z.main, names(z.main), clusterCrit::bestCriterion)
  
  # Which indices are the best with the greatest/least value?
  max.bests <- z.main %>% 
    extract(purrr::map_int(., which.max) == bests) %>% 
    cbind(z.other)
  min.bests <- z.main %>% 
    extract(purrr::map_int(., which.min) == bests) %>% 
    purrr::map_df(~ sum(.x) - .x)
  all.bests <- cbind(max.bests, min.bests)
  
  # Rank indices, sum ranks, get quantile of rank from total, filter quantiles
  ind.ranks <- purrr::map_df(all.bests, rank)
  sum.ranks <- data.frame(Algorithms = alg.all,
                          Sumrank = rowSums(ind.ranks)) %>% 
    mutate(Quantile = Sumrank / max(Sumrank)) %>% 
    filter(Quantile >= quantile)
  
  # Which algorithms are kept and which are removed? Create trimmed ensemble
  alg.keep <- as.character(sum.ranks$Algorithms)
  alg.remove <- as.character(alg.all[!(alg.all %in% alg.keep)])
  cc.trimmed <- cc.obj[, , alg.keep, k, drop = FALSE]
  
  # Reweigh only if specified and more than 1 algorithm is kept
  if (reweigh && length(alg.keep) > 1) {
    
    # Filter after knowing which to keep
    ak <- match(alg.keep, alg.all)
    max.bests <- max.bests[ak, ]
    min.bests <- z.main %>% 
      extract(ak, purrr::map_int(., which.min) == bests) %>% 
      purrr::map_df(~ sum(.x) - .x)
    
    # Create multiples of each algorithm proportion to weight
    # Divide multiples by greatest common divisor to minimize number of copies
    multiples <- cbind(as.matrix(max.bests), as.matrix(min.bests)) %>% 
      prop.table(2) %>% 
      rowMeans() %>% 
      multiply_by(100) %>% 
      round(0) %>% 
      divide_by(Reduce("gcd", .)) %>% 
      set_names(alg.keep)
    
    # Generate multiples for each algorithm, adding back dimnames metadata
    cc.trimmed <- purrr::array_branch(cc.trimmed, c(3, 4)) %>% 
      purrr::map2(., multiples, ~ rep(list(.x), .y)) %>% 
      purrr::map(abind::abind, along = 3) %>% 
      abind::abind(along = 3) %>% 
      abind::abind(along = 4) %>% 
      setter::set_dimnames(
        list(NULL,
             dimnames(cc.trimmed)[[2]],
             purrr::flatten_chr(purrr::map2(names(multiples), multiples, rep)),
             k))
  }
  # Show evaluation output
  if (show.eval) {
    eval <- z
  } else {
    eval <- NULL
  }
  return(list(alg.keep = alg.keep,
              alg.remove = alg.remove,
              eval = eval,
              data.new = cc.trimmed))
}

#' Recursively find the greater common divisor of two numbers
#' @noRd
gcd <- function(x, y) {
  r <- x %% y
  return(ifelse(r, gcd(y, r), y))
}
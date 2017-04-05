#' @param E consensus object from \code{consensus_evaluate}
#' @param n an integer specifying the top \code{n} algorithms to keep after 
#'   trimming off the poor performing ones using Rank Aggregation. If the total
#'   number of algorithms is less than \code{n} no trimming is done.
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
#'   Potentially no trimming depending on \code{n}}.
#' @rdname consensus_combine
#' @export
consensus_trim <- function(E, ii, k, n = 5, reweigh = FALSE) {
  k <- as.character(min(k))
  zk <- ii[[k]]
  alg.all <- dimnames(E)[[3]]
  
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
    extract(purrr::map_int(., which.min) == bests)
  
  # Determine trimmed ensemble using rank aggregation, only if there are more
  # algorithms than we want to keep
  if (length(alg.all) <= n) {
    alg.keep <- alg.all
  } else {
    max.ranks <- purrr::map_df(max.bests, ~rank(-.x))
    min.ranks <- purrr::map_df(min.bests, rank)
    rank.agg <- cbind(max.ranks, min.ranks) %>% 
      t() %>% 
      apply(1, function(x) alg.all[x]) %>% 
      t() %>% 
      RankAggreg::RankAggreg(., ncol(.), rho = 0.5, verbose = FALSE)
    alg.keep <- rank.agg$top.list[seq_len(n)]
  }
  alg.remove <- as.character(alg.all[!(alg.all %in% alg.keep)])
  E.trim <- E[, , alg.keep, k, drop = FALSE]
  
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
    E.trim <- purrr::array_branch(E.trim, c(3, 4)) %>% 
      purrr::map2(., multiples, ~ rep(list(.x), .y)) %>% 
      purrr::map(abind::abind, along = 3) %>% 
      abind::abind(along = 3) %>% 
      abind::abind(along = 4)
    dimnames(E.trim) <-
      list(NULL,
           dimnames(E.trim)[[2]],
           purrr::flatten_chr(purrr::map2(names(multiples), multiples, rep)),
           k)
  }
  return(list(alg.keep = alg.keep,
              alg.remove = alg.remove,
              data.new = E.trim))
}

#' Recursively find the greater common divisor of two numbers
#' @noRd
gcd <- function(x, y) {
  r <- x %% y
  return(ifelse(r, gcd(y, r), y))
}
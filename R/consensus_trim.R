#' @param quantile if an algorithm's sum of ranks across indices is in the 
#'   \code{quantile} quantile, the algorithm is kept. Otherwise it is removed 
#'   (trimmed).
#' @param reweigh logical; if \code{TRUE}, after trimming out poor performing
#'   algorithms, each algorithm is reweighed depending on its internal indices.
#' @inheritParams consensus_evaluate
#' @return \code{consensus_trim} returns a list with three elements 
#'   \item{alg.keep}{algorithms kept} \item{alg.remove}{algorithms removed} 
#'   \item{E_trimmed}{A trimmed version of a \code{ConClust} object. Potentially
#'   no different than original depending on \code{quantile}}
#' @rdname consensus_combine
#' @export
consensus_trim <- function(data, k, ..., cons.cl = NULL, ref.cl = NULL,
                           quantile = 0.75, reweigh = FALSE) {
  Sumrank <- Quantile <- NULL
  cc.obj <- abind::abind(list(...), along = 3)
  z <- consensus_evaluate(data = data, k = k, cc.obj, cons.cl = cons.cl,
                          ref.cl = ref.cl, plot = FALSE)
  alg.all <- z$internal$Algorithms
  z.main <- z$internal[, -c(1, 13:14)]
  z.other <- z$internal[, -c(1:12)]
  bests <- mapply(function(z, n) clusterCrit::bestCriterion(z, n),
                  z = as.list(z.main), n = names(z.main))
  max.bests <- cbind(z.main[, apply(z.main, 2, which.max) == bests], z.other)
  min.bests <- z.main[, apply(z.main, 2, which.min) == bests] %>% 
    apply(2, function(y) sum(y) - y)
  all.bests <- cbind(max.bests, min.bests)
  ind.ranks <- apply(all.bests, 2, rank)
  sum.ranks <- data.frame(Algorithms = z$internal$Algorithms,
                          Sumrank = apply(ind.ranks, 1, sum)) %>% 
    mutate(Quantile = Sumrank / max(Sumrank)) %>% 
    filter(Quantile >= quantile)
  alg.keep <- as.character(sum.ranks$Algorithms)
  alg.remove <- as.character(alg.all[!(alg.all %in% alg.keep)])
  cc.trimmed <- cc.obj[, , alg.keep, , drop = FALSE]
  
  if (reweigh) {
    max.bests <- max.bests %>% 
      extract(match(alg.keep, alg.all), ) 
    min.bests <- z.main[, apply(z.main, 2, which.min) == bests] %>% 
      extract(match(alg.keep, alg.all), ) %>% 
      apply(2, function(y) sum(y) - y)
    multiples <- cbind(max.bests, min.bests) %>% 
      apply(2, function(y) y / sum(y)) %>% 
      rowMeans() %>% 
      multiply_by(100) %>% 
      round(0) %>% 
      divide_by(Reduce("gcd", .)) %>% 
      set_names(alg.keep)
    cc.trimmed <- plyr::alply(cc.trimmed, 3, .dims = TRUE) %>% 
      mapply(function(d, m) rep(list(d), m), d = ., m = multiples) %>%
      lapply(abind, along = 3) %>% 
      Reduce(function(...) abind(..., along = 3), .) %>% 
      abind(along = 4)
    dimnames(cc.trimmed)[[3]] <- unname(unlist(
      mapply(rep, names(multiples), multiples)))
    dimnames(cc.trimmed)[[4]] <- z$k
  } 
  return(list(alg.keep = alg.keep,
              alg.remove = alg.remove,
              E_trimmed = cc.trimmed))
}

#' Recursively find the greater common divisor of two numbers
#' @noRd
gcd <- function(x, y) {
  r <- x %% y
  return(ifelse(r, gcd(y, r), y))
}
#' @param quantile if an algorithm's sum of ranks across indices is in the 
#'   \code{quantile} quantile, the algorithm is kept. Otherwise it is removed
#'   (trimmed).
#' @inheritParams consensus_evaluate
#' @return \code{consensus_trim} returns a list with three elements
#'   \item{alg.keep}{algorithms kept}
#'   \item{alg.remove}{algorithms removed}
#'   \item{E_trimmed}{A trimmed version of a \code{ConClust} object. Potentially
#'   no different than original depending on \code{quantile}}
#' @rdname consensus_combine
#' @export
consensus_trim <- function(data, k, ..., cons.cl = NULL, ref.cl = NULL,
                           quantile = 0.75) {
  Sumrank <- Quantile <- NULL
  cc.obj <- abind::abind(list(...), along = 3)
  z <- consensus_evaluate(data = data, k = k, cc.obj, cons.cl = cons.cl,
                          ref.cl = ref.cl, plot = FALSE)
  alg.all <- z$internal$Algorithms
  z.main <- z$internal[, -c(1, 13:14)]
  z.other <- z$internal[, -c(1:12)]
  main.names <- names(z.main)
  bests <- mapply(function(z, n) clusterCrit::bestCriterion(z, n),
                  z = as.list(z.main), n = main.names)
  max.bests <- cbind(z.main[, apply(z.main, 2, which.max) == bests], z.other)
  min.bests <- z.main[, apply(z.main, 2, which.min) == bests]
  ind.ranks <- apply(cbind(-max.bests, min.bests), 2, rank)
  sum.ranks <- data.frame(Algorithms = z$internal$Algorithms,
                          Sumrank = apply(ind.ranks, 1, sum)) %>% 
    mutate(Quantile = Sumrank / max(Sumrank)) %>% 
    filter(Quantile >= quantile)
  alg.keep <- as.character(sum.ranks$Algorithms)
  alg.remove <- as.character(alg.all[!(alg.all %in% alg.keep)])
  cc.trimmed <- cc.obj[, , alg.keep, , drop = FALSE]
  return(list(alg.keep = alg.keep,
              alg.remove = alg.remove,
              E_trimmed = cc.trimmed))
}
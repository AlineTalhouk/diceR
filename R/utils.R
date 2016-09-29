# magrittr placeholder
globalVariables(".")

#' Prepare data for consensus clustering
#'
#' Remove variables with low signal and scale before consensus clustering
#'
#' The \code{min.sd} argument is used to filter the feature space for only
#' highly variable features. Only features with a standard deviation across all
#' samples greater than \code{min.sd} will be used.
#'
#' @param data data matrix. Columns are samples and rows are genes/features.
#' @param min.sd minimum standard deviation threshold. See details.
#' @return dataset prepared for usage in \code{ConClust}
#' @author Derek Chiu
#' @importFrom stats sd
#' @export
#' @examples
#' set.seed(2)
#' x <- replicate(10, rnorm(100))
#' prepare_data(x)
prepare_data <- function(data, min.sd = 1) {
  dat.out <- data %>%
    as.data.frame() %>%
    select(which(sapply(., class) == "numeric")) %>%
    extract(apply(., 1, function(x) sd(x, na.rm = TRUE)) > min.sd,
            apply(., 2, function(x) !any(is.na(x)))) %>%
    t() %>%
    scale() %>%
    t()
  return(dat.out)
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

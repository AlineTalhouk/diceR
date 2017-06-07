#' diceR
#'
#' The diceR package performs Consensus Clustering, an ensemble clustering
#' algorithm that pools together results from different clustering realizations.
#' Each realization uses a unique combination of observations and diverse
#' clustering algorithms.
#'
#' @name diceR-package
#' @useDynLib diceR, .registration = TRUE
#' @docType package
#' @import ggplot2
#' @importFrom magrittr "%>%"
#' @importFrom purrr "%||%"
#' @importFrom Rcpp sourceCpp
NULL

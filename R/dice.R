#' Runs through the full diceR algorithm
#' 
#' A function that runs through the full diceR algorithm 
#' 
#' @param data a dataset with rows as observations, columns as variables
#' @param nk a rangeof cluster sizes (or a single value).
#' @param algorithms clustering algorithms to be used in the ensemble. Current options are c("nmfDiv", "nmfEucl", "hcAEucl", "hcDianaEucl", "kmEucl","kmSpear", "pamEucl", "pamSpear", "apEucl", "scRbf", "gmmBIC","biclust"). Please see ConClust for definition.
#' @param consensusFUNS Consensus Functions to use. Current options are majority voting, k-modes, CSPA, LCE
#' @param prune a logical value to indicate whether algorithm pruning is needed prior to consensus clustering. defaults to FALSE.
#' @param weigh a logical value to indicate whether after pruning, certain algorithms should have more weight than others. defaults to FALSE.If weigh is TRUE pruning will be set to TRUE as well.
#' @param evaluate indicate whether internal evaluation is needed. evaluate=c("internal","external","both"). If evaluate is either external or both, then a reference class must be provided.
#' @param refClass a vector of length n, indicating the reference class to compare against.

#' @export

dice <- function(data, nk,
                 algorithms = c("hcAEucl","kmEucl","scRbf", "gmmBIC"),
                 consensusFUNS = c("k-modes","CSPA","majority","LCE"),
                 prune = FALSE, weigh = FALSE, 
                 evaluate = "internal", refClass = NULL){
if(prune==FALSE & weigh==TRUE){prune <- TRUE}
if((evaluate=="external" |evaluate=="both")& refClass=NULL) {stop("Please define Reference Class or set evaluate to 'internal'")}

  }
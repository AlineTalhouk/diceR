#' Similarity to Distance
#' 
#' Converts similarity values to distance values and change matrixs format from
#' square to vector (input format for linkage function)
#' 
#' @param S N by N similarity matrix
#' @references MATLAB function stod in package LinkCluE by Simon Garrett   
#' @return distance vector
#' @author Johnson Liu
#' @export
#' 
#' @examples
#' set.seed(1)
#' E<-matrix(rep(sample(1:4,1000,replace = TRUE)),nrow=100,byrow=FALSE)
#' d <- stod(cts(E,0.8))
stod <- function(S) {
  assertthat::assert_that(is.numeric(S), is.matrix(S))
  s <- NULL
  for (a in 1:(ncol(S) - 1)) {
    s <- append(s, S[a, (a + 1):ncol(S)])
  }
  vec <- NULL
  for (i in 1:length(s)) {
    vec <- append(vec, s[[i]])
  }
  return(1 - vec)
}
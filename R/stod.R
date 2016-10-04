#' Title converts similarity values to distance values and change matrix's format from square to vector (input format for linkage function)
#'
#' @param S : N by N similarity matrix
#'
#' @return d: distance vector
#' @export
#'
#' @examples
#' data("E_LCE")
#' ASRS<-asrs(E_LCE,0.8)
#' d<-stod(ASRS)
stod<-function(S){
  assertthat::assert_that(is.numeric(S))
  assertthat::assert_that(is.matrix(S))
  s<-NULL
  for(a in 1:(ncol(S)-1)){
    s<-append(s,S[a,(a+1):ncol(S)])
  }
  vec<-NULL
  for(i in 1:length(s)){
    vec<-append(vec,s[[i]])
  }
  return(1-vec)
}

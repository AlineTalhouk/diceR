#' Title converts similarity values to distance values and change matrix's format from square to vector (input format for linkage function)
#'
#' @param S : N by N similarity matrix
#'
#' @return d: distance vector
#' @export
#'
#' @examples
stod<-function(S){
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

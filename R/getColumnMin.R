#' Title: function to find the minimum in each column of a mtrix
#'
#' @param M : input matrix
#'
#' @return vec: vector containing the minimum in each column of input matrix
#' @export
#'
#' @examples
getColumnMin<-function(M){
  vec<-NULL
  for(i in 1:ncol(M)){
    vec<-append(vec,min(M[,i]))
  }
  return(vec)
}

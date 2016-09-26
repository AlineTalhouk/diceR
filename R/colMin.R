#' Function to find the minimum in each column of a matrix
#'
#' @param M : a matrix
#'
#' @return vec: a vector containning the minimum in each column of input matrix
#' @export
#'
#' @examples
colMin<-function(M){
  vec<-NULL
  for (i in 1:ncol(M)){
    vec<-append(vec,min(M[,i]))
  }
  return(vec)
}

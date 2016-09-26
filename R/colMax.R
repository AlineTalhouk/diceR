#' Function to find the maximum in each column of a matrix
#'
#' @param M : a matrix
#'
#' @return vec: a vector containning the maximum in each column of input matrix
#' @export
#'
#' @examples
colMax<-function(M){
  vec<-NULL
  for (i in 1:ncol(M)){
    vec<-append(vec,max(M[,i]))
  }
  return(vec)
}

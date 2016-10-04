#' Function to find the minimum in each column of a matrix
#'
#' @param M : a matrix
#'
#' @return vec: a vector containning the minimum in each column of input matrix
#' @export
#'
#' @examples
#' colMin(matrix(c(1,2,3,4,5,6,7,8,9),ncol=3))
colMin<-function(M){
  assertthat::assert_that(is.numeric(M))
  assertthat::assert_that(ncol(M)>=1)
  vec<-NULL
  for (i in 1:ncol(M)){
    vec<-append(vec,min(M[,i]))
  }
  return(vec)
}

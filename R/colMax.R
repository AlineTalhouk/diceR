#' Function to find the maximum in each column of a matrix
#'
#' @param M : a matrix
#'
#' @return vec: a vector containning the maximum in each column of input matrix
#' @export
#'
#' @examples
#' colMax(matrix(c(1,2,3,4,5,6,7,8,9),ncol=3))
colMax<-function(M){
  assertthat::assert_that(ncol(M)>=1)
  assertthat::assert_that(is.numeric(M))
  vec<-NULL
  for (i in 1:ncol(M)){
    vec<-append(vec,max(M[,i]))
  }
  return(vec)
}

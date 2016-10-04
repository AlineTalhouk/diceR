#' Title function to check if a number is a positive integer
#'
#' @param x : a number
#'
#' @return TRUE if x is a positive integer and FALSE otherwise
#' @export
#'
#' @examples
checkPosInt<-function(x){
  assert_that(is.numeric(x))
  assert_that(length(x)==1)
  if(x<0){
    return(FALSE)
  } else{
    if(ceiling(x)==floor(x)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
}

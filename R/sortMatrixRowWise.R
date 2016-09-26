#' Title function to rows of a matrix
#'
#' @param M : input matrix
#' @param order : "ascending" for ascending order or "descending" for descending order
#'
#' @return a matrix whose rows are sorted accordingly
#' @export
#'
#' @examples
sortMatrixRowWise<-function(M,order){
  for(i in 1:nrow(M)){
    temp<-M[i,]
    if(order=="ascending"){
      temp<-sort(temp)
    }else if(order=="descending"){
      temp<-sort(temp,decreasing = TRUE)
    } else{
      stop("Invalid order argument in sortMatrixRowWise")
    }
    M[i,]<-temp
  }
  return(M)
}

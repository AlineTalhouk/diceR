#' Title function to get row number of the minimum in each column of a matrix
#'
#' @param M : a matrix
#'
#' @return vector containing row numbers of the minimum in each column of a matrix
#' @export
#'
#' @examples
getRowOfMin<-function(M){
  rows<-NULL
  for(i in 1:ncol(M)){
    rows<-append(rows,which(M[,i]==min(M[,i])))
  }
  return(rows)
}

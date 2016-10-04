#' Title: function to find the locations of an enetry in a matrix
#'
#' @param M : input matrix
#' @param toFind : number to be found
#'
#' @return: rows: vector containing the row numbers of toFind in M
#'          cols: vector containing the column numbers of toFind in M
#' @export
#'
#' @examples
#' getRowColNumbers(matrix(1:16,ncol=4),11)
getRowColNumbers<-function(M,toFind){
  assertthat::assert_that(is.matrix(M))
  assertthat::assert_that(is.numeric(M))
  assertthat::assert_that(is.numeric(toFind))
  assertthat::assert_that(toFind %in% M)
  assertthat::assert_that(length(toFind)==1)
  rowVec<-NULL
  colVec<-NULL
  spots<-which(M==toFind)
  for(i in 1:length(spots)){
    if(spots[i]%%nrow(M)==0){
      colVec<-append(colVec,as.integer(spots[i]/nrow(M)))
      rowVec<-append(rowVec,as.integer(nrow(M)))
    }else{
      colVec<-append(colVec,as.integer(ceiling(spots[i]/nrow(M))))
      rowVec<-append(rowVec,as.integer(spots[i]%%nrow(M)))
    }
  }
  return(list(rows=rowVec,cols=colVec))
}

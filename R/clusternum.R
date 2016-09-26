#' Title function to assign leaves under cluster c to c, equivalent to MATLAB function clusternum
#'
#' @param X
#' @param k
#' @param c
#'
#' @return
#' @export
#'
#' @examples
clusternum<-function(X,resultT,k,c){
  m<-nrow(X)+1
  children<-NULL
  while(length(k)!=0){
    children<-X[k,1:2]
    resultT[children[which(children<=m)]]<-c
    k<-children[which(children>m)]-m
  }
  return(resultT)
}

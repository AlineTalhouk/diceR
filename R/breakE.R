#' Title function to create a cluster matrix for all algorithms in ConClust during 1 of the repetitions
#'
#' @param E : ensemble out of crEnsemble using scheme=3
#' @param methods : algorithm names as those in ConClust
#' @param col_ind : indicator of which repetition. Range is 1 to dim(E)[2]
#'
#' @return a 2D matrix
#' @export
#'
#' @examples
breakE<-function(E,methods,col_ind){
  E_out<-NULL
  for(k in 1:dim(E)[3]){
    E_out<-cbind(E_out,E[,col_ind,k])
  }
  colnames(E_out)<-methods
  return(E_out)
}

#' Title Part of valid_RandIndex
#'
#' @param Mem1 : a vector
#' @param Mem2 : a vector
#'
#' @return See matlab valid_RandIndex.m in LinkCluE
#' @export
#'
#' @examples
Contingency<-function(Mem1,Mem2){
  if(!is.vector(Mem1) | !is.vector(Mem2)){
    stop("RandIndex: Requires two vector arguments")
  }
  Cont<-zeros(max(Mem1),max(Mem2))
  for(i in 1:length(Mem1)){
    Cont[Mem1[i],Mem2[i]]<-Cont[Mem1[i],Mem2[i]]+1
  }  
  return(Cont)
}

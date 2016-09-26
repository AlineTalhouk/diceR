#' Title convert one cluster result to another cluster with same algorithm but different labeling
#'
#' @param numClust : number of clusters
#' @param v2 : reference vector
#' @param vec : vector to be converted
#'
#' @return converted vec
#' @export
#'
#' @examples checkVec(4,c(4,2,1,3),c(1,3,2,4)) returns 3,4,2,1
checkVec<-function(numClust,v2,vec){
  v1<-1:numClust
  tempIndex<-0
  for(i in 1:length(vec)){
    tempIndex<-which(v2==vec[i])
    vec[i]<-v1[tempIndex]
  }
  return(vec)
}

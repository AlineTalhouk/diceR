#' Title function for computing compactness validity index for a clustering result
#'
#' @param X : a dataset, rows are observations, columns are variables
#' @param labels : cluster labels from a clustering result(a vector)
#'
#' @return compactness score
#' @export
#'
#' @examples
valid_compactness<-function(X,labels){
  n<-length(labels)
  C<-sort(unique(labels))
  k<-length(C)
  cp<-0
  for(i in 1:k){
    ind<-which(labels==C[i])
    nk<-length(ind)
    if(nk<=1){
      cp<-cp+0
    } else{
      sum_d<-0
      sum_d<-sum(dist(X[ind,],method="euclidean"))
      cp<-cp+(nk*(sum_d/(nk*(nk-1)/2)))
    }
  }
  return(cp/n)
}

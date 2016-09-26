#' Title function for computing classification accuracy for a clustering result
#'
#' @param labels 
#' @param truelabels 
#'
#' @return classification accuracy score
#' @export
#'
#' @examples
valid_CA<-function(labels,truelabels){
  nrow<-length(truelabels)
  C<-sort(unique(labels))
  k<-length(C)
  clusters<-list(NULL)
  for(i in 2:k){
    clusters<-c(clusters,NULL)
  }
  ca<-0
  for(i in 1:k){
    ind<-which(labels==C[i])
    clusters[[i]]<-matrix(ind,ncol=1)
    n<-length(ind)
    temp<-NULL
    for(j in 1:n){
      temp<-append(temp,truelabels[ind[j]])
    }
    clusters[[i]]<-cbind(clusters[[i]],temp)
    TC<-sort(unique(clusters[[i]][,2]))
    kTC<-length(TC)
    ind<-NULL
    for(l in 1:kTC){
      ind<-append(ind,length(which(clusters[[i]][,2]==TC[l])))
    }
    M<-max(ind)
    I<-which(ind==M)
    ca<-ca+M
  }
  return(ca/nrow)
}

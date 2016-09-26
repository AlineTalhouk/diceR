#' Title function to create clusters using distance criterion from a hierarchical cluster tree
#'
#' @param Z : matrix representing a hierarchical cluster tree
#' @param maxclust : desired number of clusters
#'
#' @return resultT: clustering result
#' @export
#'
#' @examples
cluster<-function(Z,maxclust){
  m<-nrow(Z)+1
  resultT<-zeros(m,length(maxclust))
  for(j in 1:length(maxclust)){
    if(m<=maxclust[j]){
      resultT[,j]<-1:m
    }else if(maxclust[j]==1){
      resultT[,j]<-rep(1,m)
    } else{
      clsnum<-1
      for(k in (m+1-maxclust[j]):(m-1)){
        i<-Z[k,1]
        if(i<=m){
          resultT[i,j]<-clsnum
          clsnum<-clsnum+1
        }else if(i<(2*m-maxclust[j]+1)){
          resultT[,j]<-LinkCluE::clusternum(Z,resultT[,j],i-m,clsnum)
          clsnum<-clsnum+1
        }
        i<-Z[k,2]
        if(i<=m){
          resultT[i,j]<-clsnum
          clsnum<-clsnum+1
        }else if(i<(2*m-maxclust[j]+1)){
          resultT[,j]<-LinkCluE::clusternum(Z,resultT[,j],i-m,clsnum)
          clsnum<-clsnum+1
        }
      }
    }
  }
  return(resultT)
}

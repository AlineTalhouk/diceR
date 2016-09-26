#' Title function for computing Davies-Bouldin index and Dunn index
#'
#' @param X : a data set whose rows are observations, columns are variables
#' @param labels : cluster labels from a clustering result(a vector)
#'
#' @return DB: Davies-Bouldin score, Dunn: Dunn score
#' @export
#'
#' @examples
valid_DbDunn<-function(X,labels){
  nrow<-nrow(X)
  nc<-ncol(X)
  k<-max(labels)
  temp<-valid_sumsqures(X,labels,k)
  st<-temp$Tot
  sw<-temp$W
  sb<-temp$B
  cintra<-temp$Sintra
  cinter<-temp$Sinter
  R<-zeros(k)
  dbs<-zeros(1,k)
  for(i in 1:k){
    for(j in (i+1):k){
      if(j<=k){
        if((cinter[i,j]==0)){
          R[i,j]<-0
        }else{
          R[i,j]<-(cintra[i]+cintra[j])/cinter[i,j]
        }
      }
    }
    dbs[1,i]<-max(R[i,])
  }
  DB<-mean(dbs[1,1:k-1])
  dbs<-max(cintra)
  R<-cinter/dbs
  for(i in 1:(k-1)){
    S<-R[i,(i+1):k]
    dbs[i]<-min(S)
  }
  return(list(DB=DB,Dunn=min(dbs)))
}

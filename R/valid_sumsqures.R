#' Title function to compute within group, between group, and total sum of squares and cross-products
#'
#' @param data : a matrix with each column representing a variable
#' @param labels : a vector indicating class labels
#' @param k : number of clusters
#'
#' @return W: within group sum of squares and cross-products; B: between group sum of squares and cross-products;T: total sum of squares and cross-products;
#'          Sintra & Sinter: centroid diameter and linkage distance
#' @export
#'
#' @examples
valid_sumsqures<-function(data,labels,k){
  data<-as.matrix(data)
  ncase<-nrow(data)
  m<-ncol(data)
  Dm<-t(colMeans(data))
  Dm<-data-Dm[ones(ncase,1),]
  Tot<-t(Dm)%*%Dm
  W<-matrix(rep(0,nrow(Tot)*ncol(Tot)),nrow=nrow(Tot))
  Dm<-matrix(rep(0,k*m),nrow = k)
  Sintra<-matrix(rep(0,k),nrow=1)
  for(i in 1:k){
    if(k>1){
      Cindex<-which(labels==i)
    } else{
      Cindex<-1:ncase
    }
    nk<-length(Cindex)
    if(nk>1){
      dataC<-data[Cindex,]
      m<-colMeans(dataC)
      Dm[i,]<-m
      dataC<-dataC-repmat(m,nk,1)
      W<-W+t(dataC)%*%dataC
      dataC<-rowSums(dataC^2)
      Sintra[i]<-mean(sqrt(dataC))
    }
  }
  B<-Tot-W
  Sinter<-matrix(rep(0,k^2),nrow=k)
  if(k>1){
    for(i in 1:k){
      for(j in (i+1):k){
        if(j<=k){
          m<-abs(Dm[i,]-Dm[j,])
          Sinter[i,j]<-sqrt(sum(m^2))
          Sinter[j,i]<-Sinter[i,j]
        }
      }
    }
  }
  return(list(Tot=Tot,W=W,B=B,Sintra=Sintra,Sinter=Sinter))
}

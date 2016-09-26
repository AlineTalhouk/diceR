#' Title A function for computing approximated simrank based similarity matrix
#'
#' @param E : matrix of cluster ensemble
#' @param dc :decay factor, range 0 to 1, inclusive
#'
#' @return S: ASRS matrix
#' @export
#'
#' @examples
asrs<-function(E,dc){
  n<-nrow(E)
  M<-ncol(E)
  E.new<-relabelCl(E)
  E<-E.new$newE
  no_allcl<-E.new$no_allcl
  wcl<-weightCl(E)
  CS<-matrix(rep(0,no_allcl*no_allcl),nrow=no_allcl)
  for(i in 1:(no_allcl-1)){
    Ni<-wcl[i,]
    ni<-length(Ni[which(Ni>0)])
    for(j in (i+1):no_allcl){
      Nj<-wcl[j,]
      nj<-length(Nj[which(Nj>0)])
      if((ni*nj)>0){
        CS[i,j]<-(Ni%*%Nj)/(ni*nj)
      }
    }
  }
  if(max(CS)>0){
    CS<-CS/max(CS)
  }
  CS<-CS+t(CS)
  for(i in 1:no_allcl){
    CS[i,i]<-1
  }
  S<-matrix(rep(0,n*n),nrow=n)
  for(i in 1:(n-1)){
    for(ii in (i+1):n){
      for(j in 1:M){
        for(jj in 1:M){
          if(CS[E[i,j],E[ii,jj]]==1){
            S[i,ii]<-S[i,ii]+1
          } else{
            S[i,ii]<-S[i,ii]+dc*CS[E[i,j],E[ii,jj]]
          }
        }
      }
    }
  }
  S<-S/(M*M)
  S<-S+t(S)
  for(i in 1:n){
    S[i,i]<-1
  }
  return(S)
}

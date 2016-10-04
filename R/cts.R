#' Title : fucntion for computing connected triple based similarity matrix
#'
#' @param E : matrix of cluster ensemble
#' @param dc : decay factor, [0, 1]
#'
#' @return S: CTS matrix
#' @export
#'
#' @examples 
#' data("E_LCE")
#' CTS<-cts(E_LCE,0.8)
cts<-function(E,dc){
  assertthat::assert_that(is.matrix(E))
  assertthat::assert_that(dc>=0 && dc<=1)
  n=nrow(E)
  M=ncol(E)
  E.new<-relabelCl(E)
  E<-E.new$newE
  no_allcl<-E.new$no_allcl
  wcl<-weightCl(E)

  wCT<-matrix(rep(0,no_allcl*no_allcl),nrow=no_allcl)

  maxCl<-NULL
  minCl<-NULL

  for(i in 1:ncol(E)){
    maxCl<-append(maxCl,max(E[,i]))
    minCl<-append(minCl,min(E[,i]))
  }

  for(q in 1:M){
    for(i in minCl[q]:(maxCl[q]-1)){
      Ni<-wcl[i,]
      for(j in (i+1):(maxCl[q])){
        Nj<-wcl[j,]
        wCT[i,j]<-sum(colMin(rbind(Ni,Nj)))
      }
    }
  }
  if(max(max(wCT))>0){
    wCT<-wCT/max(max(wCT))
  }
  wCT<-wCT+t(wCT)
  for(i in 1:no_allcl){
    wCT[i,i]<-1
  }
  S<-matrix(rep(0,n*n),nrow=n)
  for(m in 1:M){
    for(i in 1:(n-1)){
      for(j in (i+1):n ){
        if(E[i,m]==E[j,m]){
          S[i,j]=S[i,j]+1
        } else{
          S[i,j]=S[i,j]+dc*wCT[E[i,m],E[j,m]]
        }
      }
    }
  }
  S<-S/M
  S<-S+t(S)
  for(i in 1:n){
    S[i,i]<-1
  }
  return(S)
}

#' Title function for computing simrank based similarity matrix
#'
#' @param E : N by M matrix of cluster ensemble
#' @param dc : decay factor, [0,1]
#' @param R : number of iterations for simrank algorithm
#'
#' @return N by N SRS matrix
#' @export
#'
#' @examples
#' data("E_LCE")
#' SRS<-srs(LCE,0.8,5)
srs<-function(E,dc,R){
  assertthat::assert_that(is.matrix(E))
  assertthat::assert_that(dc>=0 && dc<=1)
  assertthat::assert_that(is.numeric(R))
  assertthat::assert_that(checkPosInt(R)==TRUE)
  n<-nrow(E)
  M<-ncol(E)
  E.new<-relabelCl(E)
  E<-E.new$newE
  no_allcl<-E.new$no_allcl
  S<-diag(x=1,nrow=n,ncol=n)
  C<-diag(x=1,nrow=no_allcl,ncol=no_allcl)

  for(r in 1:(R-1)){
    S1<-diag(x=1,nrow=n,ncol=n)
    for(i in 1:(n-1)){
      Ni<-E[i,]
      for(ii in (i+1):n){
        sum_sim<-0
        Nii<-E[ii,]
        for(k in 1:M){
          for(kk in 1:M){
            sum_sim<-sum_sim+C[Ni[k],Nii[kk]]
          }
        }
        S1[i,ii]<-(dc/(M*M))*sum_sim
      }
    }
    S1<-S1+t(S1)
    for(i in 1:n){
      S1[i,i]<-1
    }
    C1<-diag(x=1,nrow=no_allcl,ncol=no_allcl)
    for(i in 1:(no_allcl-1)){
      Ni<-getRowColNumbers(E,i)$rows
      col<-getRowColNumbers(E,i)$cols
      nki<-length(Ni)
      for(ii in (i+1):no_allcl){
        sum_sim<-0
        Nii<-getRowColNumbers(E,ii)$rows
        col<-getRowColNumbers(E,ii)$cols
        nkii<-length(Nii)
        for(k in 1:nki){
          for(kk in 1:nkii){
            sum_sim<-sum_sim+S[Ni[k],Nii[kk]]
          }
        }
        if((nki*nkii)>0){
          C1[i,ii]<-dc*sum_sim/(nki*nkii)
        }
      }
    }
    C1<-C1+t(C1)
    for(i in 1:no_allcl){
      C1[i,i]<-1
    }
    S<-S1
    C<-C1
  }
  return(S)
}

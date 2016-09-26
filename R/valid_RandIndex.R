#' Title function to calculate Rand indices to compare two partitions
#'
#' @param c1 : vector listing class membership
#' @param c2 : another vector listing class membership
#'
#' @return AR: adjusted Rand index, RI: unadjusted Rand index, MI: Mirkin's index, HI: Hubert's index
#' @export
#'
#' @examples
valid_RandIndex<-function(c1,c2){
  if(!is.vector(c1) | !is.vector(c2)){
    stop("RandIndex: Requires two vector arguments")
  }
  C<-Contingency(c1,c2)
  n<-sum(C)
  nis<-sum(rowSums(C)^2)
  njs<-sum(colSums(C)^2)
  t1<-n*(n-1)/2
  t2<-sum(C^2)
  t3<-0.5*(nis+njs)
  nc<-(n*(n^2+1)-(n+1)*nis-(n+1)*njs+2*(nis*njs)/n)/(2*(n-1))
  A<-t1+t2-t3
  D<-t3-t2
  if(t1==nc){
    AR<-0
  }else{
    AR<-(A-nc)/(t1-nc)
  }
  RI<-A/t1
  MI<-D/t1
  HI<-(A-D)/t1
  return(list(AR=AR,RI=RI,MI=MI,HI=HI))
}

#' Title function to create hierarchical cluster tree
#'
#' @param Y : a distance vector
#' @param method : linkage method, can be "average","single", or "complete"
#'
#' @return matrix denoting hierarchical cluster tree
#' @export
#'
#' @examples
linkage<-function(Y,method){
  n<-length(Y)
  m<-ceiling(sqrt(2*n))
  Z<-zeros(m-1,3)
  N<-rep(0,2*m-1)
  N[1:m]<-1
  n<-m
  R<-1:n
  if(any(c("ce","me","wa")==method)){
    Y<-Y*Y
  }
  for(s in 1:(n-1)){
    if(method=="average"){
      if((m-1)>=2){
        p<-seq(m-1,2,-1)
      }
      I<-rep(0,0.5*m*(m-1))
      I[cumsum(c(1,p))]<-1
      I<-cumsum(I)
      J<-rep(1,0.5*m*(m-1))
      J[cumsum(p)+1]<-2-p
      J[1]<-2
      J<-cumsum(J)
      W<-N[R[I]]*N[R[J]]
      W<-W[!is.na(W)]
      tempVec<-Y/W
      v<-min(tempVec)
      k<-which(tempVec==min(tempVec))[1]
    }else{
      v<-min(Y)
      k<-which(Y==min(Y))[1]
    }
    if((length(m)!=0)&(length(k)!=0)){
      i<-floor(m+1/2-sqrt(m^2-m+1/4-2*(k-1)))
      j<-k-(i-1)*(m-i/2)+i
    }
    Z[s,]<-cbind(R[i],R[j],v)
    tempUIJ<-UIJ(i,j,m)
    U<-tempUIJ$U
    I<-tempUIJ$I
    J<-tempUIJ$J

    if(!is.null(I)){
      if(method=="single"){
        Y[I]<-colMin(rbind(Y[I],Y[J]))
      } else if(method=="complete"){
        Y[I]<-colMax(rbind(Y[I],Y[J]))
      } else if(method=="average"){
        Y[I]<-Y[I]+Y[J]
      } else{
        stop("Invalid method for function linkage")
      }
    }
    J<-c(J,i*(m-(i+1)/2)-m+j)
    Y<-Y[-J]
    m<-m-1
    N[n+s]<-N[R[i]]+N[R[j]]
    R[i]<-n+s
    if(((n-1)>=j)&(n>=(j+1))){
      R[j:(n-1)]<-R[(j+1):n]
    }
  }
  Z[,1:2]<-sortMatrixRowWise(Z[,1:2],"ascending")
  return(Z)
}
